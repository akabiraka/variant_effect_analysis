import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
home_dir = "../"

import time
import os
from Bio import Entrez
from urllib.error import HTTPError
import xml.etree.ElementTree as ET
import pandas as pd
import xmltodict
from typing import List


Entrez.email = "akabir0101@gmail.com" # provide your user email 
# RECOMMENDED: apply for API key from NCBI (https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/). 
# 10 queries per second with a valid API key, otherwise 3 queries per seconds are allowed for 'None'
Entrez.api_key = "328570309ccd040632796143ec88b51bcf08"
retmax = 100 # return 20 rs per batch example, max=1000
namespace = "https://www.ncbi.nlm.nih.gov/SNP/docsum"
ns = u'{%s}' % namespace
nsl = len(ns)

def search(snp_ids: str):
    #snp_ids: '121913562 121913566 80358221',
    # dbSNP supported query terms (https://www.ncbi.nlm.nih.gov/snp/docs/entrez_help/) can be build and test online using web query builder (https://www.ncbi.nlm.nih.gov/snp/advanced) 
    # esearch handle
    eShandle = Entrez.esearch(db="snp",  # search dbSNP
                            term=snp_ids, #'121913562 121913566 80358221',
                            usehistory="y", #cache result on server for download in batches
                            retmax=retmax
                            )                        

    # get esearch result
    global eSresult, webenv, query_key, total_count
    eSresult = Entrez.read(eShandle)
    webenv = eSresult["WebEnv"]
    query_key = eSresult["QueryKey"]
    total_count = int(eSresult["Count"])
    print(f"Query result count: {total_count}, Fetch count: {len(range(0, total_count, retmax))}")


def parse_response_xml(x:str):
    o = xmltodict.parse(x) #returns dict obj
    docs = o["ExchangeSet"]["DocumentSummary"]

    outs = []
    for i, doc in enumerate(docs):
        genes = doc["GENES"]["GENE_E"] # can have multiples
        if isinstance(genes, list):
            genes = [f'{gene["NAME"]}:{gene["GENE_ID"]}' for gene in genes]
            genes = ",".join(genes)
        elif isinstance(genes, dict):
            genes = f'{genes["NAME"]}:{genes["GENE_ID"]}'
        # print(genes)

        
        mafs = doc["GLOBAL_MAFS"]["MAF"] # can have multiples
        if isinstance(mafs, list):
            mafs = [f'{maf["STUDY"]}:{maf["FREQ"]}' for maf in mafs]
            mafs = ",".join(mafs)
        elif isinstance(mafs, dict):
            mafs = f'{mafs["STUDY"]}:{mafs["FREQ"]}'
            # print(mafs)

        variations = doc["DOCSUM"] # an example: HGVS=NC_000023.11:g.154360389T>A,NC_000023.11:g.154360389T>C,NW_003871103.3:g.1794368T>A,NW_003871103.3:g.1794368T>C,NG_011506.2:g.19250A>T,NG_011506.2:g.19250A>G,NM_001456.4:c.3406A>T,NM_001456.4:c.3406A>G,NM_001456.3:c.3406A>T,NM_001456.3:c.3406A>G,NM_001110556.2:c.3406A>T,NM_001110556.2:c.3406A>G,NM_001110556.1:c.3406A>T,NM_001110556.1:c.3406A>G,NC_000023.10:g.153588757T>A,NC_000023.10:g.153588757T>C,NP_001447.2:p.Ile1136Phe,NP_001447.2:p.Ile1136Val,NP_001104026.1:p.Ile1136Phe,NP_001104026.1:p.Ile1136Val|SEQ=[T/A/C]|LEN=1|GENE=FLNA:2316
        variations = (variations[5:].split("|")[0]).split(",")
        variations = [v for v in variations if v.startswith("NP_")]
        variations = ",".join(variations)
        
        data = {
            "snp_id": doc["SNP_ID"],
            "acc": doc["ACC"],
            "chrpos": doc["CHRPOS"],
            "spdi": doc["SPDI"],
            "tax_id": doc["TAX_ID"],
            "snp_class": doc["SNP_CLASS"],
            "create_date": doc["CREATEDATE"],
            "update_date": doc["UPDATEDATE"],
            "clinical_significance": doc["CLINICAL_SIGNIFICANCE"],
            "fxn_class": doc["FXN_CLASS"],
            "validated": doc["VALIDATED"],
            "genes": genes,
            "mafs": mafs,
            "variations": variations
        }
        outs.append(data)
        
        # if i==5: break
    return outs

def save(data, out_filepath:str):
    # data: list of objects
    out_df = pd.DataFrame(data)
    if not os.path.exists(out_filepath):
        out_df.to_csv(out_filepath, chunksize=10000, sep="\t", index=False, mode="a", header=True)
    else:
        out_df.to_csv(out_filepath, chunksize=10000, sep="\t", index=False, mode="a", header=False)
        
        

def download_batch(start, retmax):
    global eSresult, webenv, query_key, total_count
    attempt = 0
    while (attempt < 3):
        attempt += 1
        try:
            fetch_handle = Entrez.efetch(db="snp",
                                        # rettype="",
                                        retmode="xml",
                                        retstart=start,
                                        retmax=retmax,
                                        webenv=webenv,
                                        query_key=query_key)
            # break
        except HTTPError as err:
            if 400 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(10)
            else:
                raise
    try:                
        data = fetch_handle.read().decode()
        data = parse_response_xml(data)
        fetch_handle.close()
        return data
    except:
        print(f"Error. Downloading again record {start+1} to {start+retmax}")
        return download_batch(start, retmax)
    # return fetch_handle
    
    
    
def download_and_save():
    # sample codes adopted with modifications from http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc139.
    out_filepath = "data/SNPdbe/snps.tsv"
    # if os.path.exists(out_filepath): os.remove(out_filepath)

    fetch_count = 1 # 1-indexed
    start_idx = 0 # 0-indexed
    for start in range(start_idx, total_count, retmax):
        end = min(total_count, start+retmax)
        print(f"Fetch no: {fetch_count}\tDownloading record {start+1} to {end}")
            
        data = download_batch(start, retmax)
        save(data, out_filepath)
            
        # if fetch_count==2: break
        fetch_count += 1            

# search('121913562 121913566 80358221')
snps_txt = '121913562 121913566 80358221 56053615 61751374 61750146 61751402 1800555 61750646 1800548 1800553 61750639 62645940 76157638 121909203 61753034 61751392 1801466 1800552 61748536 61750152 61749450 61750135 179363882 121917869 121917867 114025668 61732239 121908737 28941471 119450941 121908522 121908524 121908529 121908523 121908525 4426527 34116584 671 145078268 1800546 77718928 78340951 118204430 151052374 4547 145467699 138840536 121918008 34605986 63751122 121909571 121909556 139392083 121912713 1130409 33956927 121912727 147210663 1801689 1801690 121912442 121912438 121912432 121912431 121912436 121912434 139629762 71581996 145669559 104894897 104894888 59962885 121913000 121913002 61726464 121913001 61726465 121913003 121434566 79184941 78311289 121913105 28928868 1138272 1695 28936396 121909308 148640446 80356482 1800562 113993946 143517122 137852783 1801275 121913512 148698650 118204082 118204067 104893747 56268439 28934905 28934907 28934904 28935468 28934906 121913560 13447332 63750376 63750129 63751438 63750570 74315448 74315447 2234916 137853096 28937597 132630303 118204443 118204438 118204435 118204442 118204437 118204441 118204444 118204446 118204436 121964866 72552734 104893977 151344537 190549838 104893751 68026851 72554308 111033244 111033308 80338848 28939086 121909218 137852697 137852696 62516092 62507279 62516101 28939673 104894173 62619919 28939674 104894174 28939671 104894178 63751141 63751037 75076352 74799832 79658334 104894230 104893981 751141 121909301 121909303 72554662 63750756 63750424 63751273 137853301 137853300 137853298 121918000 121918002 121918001 4986790 4986791 137853247 121964858 121964855 121964857 121964856 121434597 113578517 104894758 104894759 121909653 121909651 121909652 55667289 121909650 104894749 104894748 120074187 74462309 28940572 74315453 28940571 104893878 104893877 74315439 28939374 121908415 28939376 121908417 28939375 121908416 104894848 121913634 121913632 4994 121434582 121909675 28934575 28934271 121912659 28934578 121912651 28934576 35414700 74315507 121907913 121918397 111033593 62638191 62638182 104894374 62638193 62638185 28934882 119103215 28934881 119103213 28934883 119103219 119103221 119103220 150591260 61749418 1800551 1762111 28938473 1800550 1800549 61754030 139655975 121912703 121909504 861539 79281338 75002628 63750526 63751235 63750306 137853104 28931573 28941472 145138923 121434255 28897696 80357064 1799950 28934895 79761867 137852872 137852870 137852873 128620183 115129687 3732880 121907950 116840789 121909276 116840773 121909277 116840805 116840778 121909278 28939087 10509681 66501115 11572080 4987161 28371759 4986910 4986913 104893706 104893701 104893708 80357410 80357498 28897672 80357406 80357000 80357163 80356929 80357382 80356880 80357017 17879961 137853007 121912443 80265967 121912433 5742905 121964970 121964971 28935180 104894896 104894899 28935481 28935482 121434622 121965074 80338897 121965073 121965078 11555096 32817201 121909672 28933070 121909673 121917877 80359814 28371560 4988496 139763309 4988498 121918120 121918119 121918118 111033780 111033715 111033773 111033687 2070074 111033721 75391579 111033735 111033754 111033800 111033658 111033690 111033817 111033701 111033796 56275071 121909737 2230288 121908311 421016 76763715 76910485 80356771 104893914 104893913 104893909 121909750 104893967 28928905 28936703 119481080 50492298 104893858 121918101 137852510 121918352 35761929 121918351 137852624 147181709 74315391 28939684 104894575 104894578 121912507 74315445 120074180 120074195 118192205 61027685 11554495 57749775 121909141 121909143 121909142 148499544 104893843 104893842 104893848 121912953 77931234 121434281 121434283 121434280 2236700 1044498 74315492 74315493 1801133 1801131 28934897 104894363 104894368 121913658 28933099 104894369 74315494 137854451 28929493 137854448 137854450 28929494 137854449 104894609 104894608 104894603 121909249 104893758 121434598 137852695 148412181 62652721 118203925 62652693 74603784 62652698 75193786 62652700 76394784 104894179 28938169 28939672 104893761 5938 28936379 63750215 63751229 63750231 661 63751309 63750543 63750301 63750004 63751416 63749962 74315403 74315411 142771326 1805076 75996173 137853259 113403872 76087194 76534745 75075748 77316810 144029481 122454125 121908572 121908571 121908573 121908574 121913295 118192170 63749869 118192150 118192176 111033627 79389353 121908484 121908479 77301881 76163360 75660264 104893926 76871093 104893925 104893922 104894160 104894158 104894161 137852677 2071203 104894778 104894771 104894772 104894768 104894769 104894773 104894777 104894381 104894378 121918687 121918697 137853130 137853129 128620184 1799899 121918079 113994174 187322901 28941775 5030824 5030804 5030821 5030827 104893824 5030809 104893826 104894756 104894760 121909798 121909796 121913400 104894504 121918447 121913507 151344462 28934571 121913343 28934573 121912662 121912664 121908842 33974228 118204095 118204101 118204096 119481075 119481077 121907987 80358220 121912744 28929477 104894338 28931580 149659001 62542743 139297434 4008659 104894328 121913448 121913449 121913459 121913461 28937590 80357365 80357013 80357353 80357492 66468541 36209567 63750959 63751165 114459867 63750416 63750635 79738788 121918394 137852238 121909607 137852482 80356668 62623459 121913529 121913530 121913254 121913233 137852494 137852485 121909547 28933979 121912719 118203931 118203932 121913496 11554290 104894229 17851045 121917756 104894226 104894228 72552709 121918102 121913250 149800596 151344534 121909551 121913088 28942085 121918477 121912712 1041981 142541719 104894453 76285851 121909550 121909554 72648363 72552710 121918478 137852230 121909606 121434595 33918343 121918393 121909555 72653134 28999113 28930970 121908714 121909557 28929468 121909558 121909552 121909553 429358 769455 7412 11548791 121913087 121913094 121913092 121909617 121909616 137852481 137852528 137852529 137852530 28931569 121913237 121918381 121909790 121909791 137852266 76901040 79527524 121909549 121964925 28940870 121907948 28930978 121909548 121912717 137852243 137852225 137852233 137852249 137852241 137852260 137852234 137852236 137852232 137852259 137852268 28935499 5907 137852479 121913091 121918480 132630279 1804495 137852264 121909605 104894831 121909533 121909563 118204057 121908739 121908738 121908740 121908736 28929469 121918392 137852265 137852240 137852242 121964939 118204056 28931568 121918479 121918445 77645174 121918395 28929476 61753185 138310841 41295257 148474991 118204061 118204062 5030869 137852322 137852315 137852313 137852314 137852323 137852317 137852320 1050829 137852329 74575103 137852321 137852316 137852324 72554665 72554664 5030868 137852326 137852327 1048994 3189053 1048991 71421642 1130341 1048992 1052753 118204063 121965006 62652713 5030851 180710258 61754393 104894830 138567132 137852274 121909610 41449150 121908050 118204064 28934893 118204059 118204060 62508588 62514927 143312735 104893908 1800458 137852275 137852433 28933673 137852458 137852397 28933679 137852395 28933672 28933675 137852472 137852422 28937298 137852453 137852461 28936970 137852464 28937277 137852416 28937278 137852417 137852460 28937269 137852406 28937275 137852412 116555717 104894135 28940585 121964926 137852254 137852269 121908048 121908049 28933688 1071748 17433856 17433877 1136757 113281514 10805890 104894506 5036 78574148 79047363 79228041 77544362 80002911 104894137 121909620 121909609 121909618 75848804 111033663 121918444 28934896 121917789 121965007 121965009 52795188 1050086 28933689 121918072 121918069 121918070 121964846 104894974 104894959 104894964 104894969 104894957 141614092 104893768 28940880 62645911 121908011 62645921 104894313 17415673 62645918 1126809 61754388 63749964 63750445 63751263 121965043 1800456 121965044 121965039 121965040 121965042 104893755 137852532 14926 37312035 121909567 74315293 137852376 137852375 5331 121909566 121909564 72656324 66721653 72653143 72653154 121907947 121907949 121918446 17401309 34716432 76411309 1136625 118204075 118204068 118204069 112909524 121908028 121912983 121912982 137852299 137854543 79377490 63750264 28933992 104893779 28933994 104893773 104893780 28933993 104893781 29001637 28933395 104893795 28933390 8154453 1057910 1799853 28371674 1057909 104894656 63750671 63750579 104893789 121912759 121908722 121908721 121908719 111033792 118204072 118204073 118204071 118204076 121913547 137852541 137852540 121912439 121912435 4882 121912722 121912720 121912721 78310315 121434534 118203933 72653172 72654797 72656306 72656317 72656330 72653137 150525536 121964928 28933677 28937305 137852468 28933681 137852387 28937289 28933670 28937293 137852445 28933674 28935213 137852393 28935202 111033613 28937294 137852426 137852427 28936969 137852414 28935216 137852403 137852428 137852459 28937281 137852420 28937272 28935201 137852358 28937268 137852404 28937282 137852430 137852431 104893863 104893864 121909569 17355168 74315294 2303790 33937393 35209591 33966487 34831026 33918778 79908535 121913631 121913624 78869449 137854541 137854542 28934603 137854544 121918071 28933980 76992529 72554339 72556293 28928878 104894658 121434278 28940574 28940279 4482466 60822373 28371588 180965157 121964941 137852300 121918408 108640286 72555399 1064588 74446775 1059455 74452350 41559117 79993619 74971492 1071742 41552219 2075684 80321556 1136683 75646802 1059460 76523677 16822853 268 118204078 137854550 34377097 2234922 17417482 1051740 28935478 121908097 121908098 151219882 121912444 1799895 121912437 121912441 121434529 121434530 28935496 104894757 121913016 143370956 1303 709932 28929470 6647 121918148 121918149 1042714 1042713 38123790 28936082 104894838 28935194 104894837 28935493 121907982 137852301 121909546 147676453 104894329 104894326 149707394 104893770 29001653 104893769 104893790 62625011 72551353 72551342 72551345 72551343 28939068 700519 121913020 121913026 149541409 118204039 118204037 3894326 28362459 17855739 778805 104894438 104894433 143385179 111700736 137852504 113954997 104894274 104894273 121912445 118204002 118204001 118204004 118204003 137852224 104894968 144172724 148450964 709055 1131275 1131215 707912 151341295 41562914 41551516 3180380 142318374 41562613 1131201 138659308 151341188 41548113 12697943 145937432 41540317 41542423 41563312 141484466 151341168 151341195 41562013 151341218 9266150 137966817 2596492 148606135 121918057 121918692 104893948 28941474 5030857 121918324 121918010 121918016 140549609 104893942 121964938 17683430 33954001 17683011 121917887 138742870 121918060 121918061 28935484 74462743 1141814 74731340 104894177 121909309 137853249 121912455 121912460 121913548 67120076 181030365 121918686 121918706 121918707 121918693 121913046 137852279 121909221 66724222 72558454 119481078 2066479 119481079 121964997 137852305 699 121912749 121912835 121912842 121912843 9332739 11547328 61754278 104894139 121964932 121909365 28936702 104893837 28933074 149089920 72558409 121434599 62642932 28934900 118203921 62642930 62642931 5030847 62644503 62514951 104893762 63749805 121917733 28936676 121913014 1799958 121908004 1800556 28940872 57443665 121908096 121918009 118204017 121907985 1052133 121908727 121908715 121908726 121908732 121908723 121908730 151344535 121908725 121908734 28930969 121908717 121908733 28930971 121908735 121908716 121908731 63750646 28929489 104894201 35486059 80358219 35887327 28936073 104894800 28935174 104894657 11552822 121918399 104894330 104894332 104893989 104893993 104894673 74315304 74315302 28934891 1800566 104894159 3026906 137852231 72547567 72547568 72547576 72547563 72547569 67939114 72547572 72547559 72547560 72547571 72547566 72547573 72547562 72547556 72547570 121909608 449856 104894443 104894445 121917817 121917818 41378349 33961444 34933751 33958739 33926796 137852785 137852786 137852784 1801278 121912505 121913517 104893839 104893838 104893840 1805009 11547464 1805007 1805008 1805005 61749755 121434297 121434296 104895317 104895295 104895304 104895319 104894424 121909180 72552735 121434600 63749836 63751420 63751139 63750800 63750590 63749885 63750322 63751163 63750522 63751106 63749967 63750601 72552272 121918654 34833812 121918007 121918012 121918011 138690664 3200254 63750512 121913036 121913037 121913038 17848368 2229707 1801253 113994175 113994172 120074184 120074185 120074186 80357462 28897683 41293463'

snps_list = snps_txt.split(" ")

for i in range(0, len(snps_list), 100):
    query = " ".join(snps_list[i: i+100])
    print(i, i+100)

    search(query)
    download_and_save()
    break