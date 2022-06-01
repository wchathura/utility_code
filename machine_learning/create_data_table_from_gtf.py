import gffutils
import pandas as pd

## this set of functions need  gffutils
# conda install -c bioconda gffutils
##  Read feature database
## read the feature input databse file, read and return a gffutils


## create the Sp database -> Create a database that can be read by gffutils
#need to put the absolute path for the GTF file
#fn = gffutils.example_filename('/home/cwijes1/projects/resources/sp_resources/SpV88.3_gene.sorted.gtf')
#db = gffutils.create_db(fn, dbfn='spv88_3.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)


def read_featureDB(feature_db_file):
    ''' Read a feature db file created from gffutils'''
    import gffutils
    feature_db= gffutils.FeatureDB(feature_db_file,keep_order=True)
    return feature_db

def get_transcripts_from_featureDb(feature_db):
    '''Given a gffutils feature db file read all the transcript ids and return them as a list'''
    import gffutils
    
    trans_list=[]
    for feature in feature_db.all_features(): 
        try:
            trans=feature.attributes["transcript_id"]
            trans_list.append(trans)
        except :
           # print(feature.attributes)
            continue
    set_tr_list=set([tr[0] for tr in trans_list])
    list_of_transcripts=list(set_tr_list)
    return list_of_transcripts

def extract_coordianteFor_spliceAi(trnscript_list,feature_db,file_name_to_write,feature_list):
    ## we need to pass a
    ## trnscript_list -> list of transcripts 
    ## feature_db -> gffutils feature db file read using gffutils.FeatureDB function 
    ## file_name_to_write -> file name of the output  
    ## feature_list -> which feature type to extract passed as a list 
    import gffutils
    
    list_to_write=[]
    i=0
    for feture_type in feature_list:
       # print(feture_type)
        for trans in trnscript_list:
            final_list=[]
            exon_start_list=[]
            exon_end_list=[]
            chr_id=[]
            strand=[]
            exon_number=[]


            for f in feature_db.children(trans, featuretype=feture_type, order_by='start'):
               # if "Un" not in f[0]:
                    #print(f)
                    try:
                        #print(f)
                        chr_id.append(f[0].lower())#.split("-")[0])
                        #print(chr_id)
                        strand.append(f[6])
                        exon_start_list.append(f[3])
                        exon_end_list.append(f[4])
                    except Exception as e:
                        print(e)
                        continue
           # print(exon_end_list)

            try :
                    exon_number.append(len(exon_start_list))
                    final_list.append(chr_id[0])
                    final_list.append(strand[0])
                    final_list.append(trans)
                    final_list.append(exon_number) ## get the number of genes by getting the length of the exon coor list
                    final_list.append(exon_start_list[0]) ## add the start coor of the exon as the gene start
                    final_list.append(exon_end_list[-1])  ## add the end coor of the exon as the gene end
                    final_list.append(exon_start_list)
                    final_list.append(exon_end_list)
                    list_to_write.append(final_list)
                    i=i+1
            except:
                 continue
            #print(list_to_write)
        
    import re
    csv_to_write=[]
    for lin in list_to_write:
        xx="\t".join(str(ll) for ll in lin)
        # print(re.findall(r'\[+', xx))
        xx=(re.sub(r'\[+',"", xx))
        xx=(re.sub(r'\]+',"", xx))
    
        csv_to_write.append(xx+"\n")
        
    with open(file_name_to_write+".csv","w") as at_splice_pre:
            for csv in csv_to_write:
                at_splice_pre.write(csv)