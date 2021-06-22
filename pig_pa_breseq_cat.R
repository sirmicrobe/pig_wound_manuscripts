
###### PA14 set ######
# Set working directory to input and save files
setwd("~/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/wt_set/pa14_subtractref/")

# Save all filenames in directory ending in .txt as a variable temp
temp = list.files(pattern="*.txt")
all_snp <- lapply(temp,read.delim) #create a list with each object in list a data frame of mutation table
# headers we want to keep
headers <- c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","gene_name","gene_position","gene_product","mutation_category","new_seq","position","seq_id","size","snp_type","type","title")
#take a look at the first data frame in the list
df_snp_01 <- all_snp[[1]]
df_snp_01 <- df_snp_01[,headers] #select only headers we want to keep
head(df_snp_01)
# convert all data frames in the list to one large data frame
all_pa14_snp_df <- as.data.frame(data.table::rbindlist(all_snp,fill=T),colClasses = c("character"))
all_pa14_snp_df <- all_pa14_snp_df[,headers] #keep headers we want
all_pa14_snp_df$mutation <- paste(all_pa14_snp_df$aa_ref_seq,all_pa14_snp_df$aa_position,all_pa14_snp_df$aa_new_seq, sep=":") #create mutation column
snp_pa14_df_30 <- subset(all_pa14_snp_df, title == "sample_30") #subset just sample 30
all_pa14_snp_df <- subset(all_pa14_snp_df, title != "sample_30") #get rid of sample that makes data frame too large (mutator or off reference)
all_pa14_snp_df$strain <- "PA14"
View(all_pa14_snp_df)
head(snp_pa14_df_30)

###### PA01 set #######

setwd("~/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/wt_set/pa01_subtractref/")
temp_pa01 = list.files(pattern="*.txt")
temp_pa01 = c( "03_pa01_subtractref.txt","12_pa01_subtractref.txt","19_pa01_subtractref.txt","25_pa01_subtractref.txt", "26_pa01_subtractref.txt")
all_pa01_snp <- lapply(temp_pa01,read.delim)
all_pa01_snp
all_pa01_snp_df <- as.data.frame(data.table::rbindlist(all_pa01_snp,fill=T),stringsAsFactors=F,colClasses = c("character"))
all_pa01_snp_df <- all_pa01_snp_df[,headers] #keep headers we want
all_pa01_snp_df$mutation <- paste(all_pa01_snp_df$aa_ref_seq,all_pa01_snp_df$aa_position,all_pa01_snp_df$aa_new_seq, sep=":") #create mutation column
all_pa01_snp_df$strain <- "PA01"
View(all_pa01_snp_df)
#add pa01 calls to end of pa14 calls data frame 
#all_snp_df <- rbind(all_pa14_snp_df, all_pa01_snp_df)
#View(all_snp_df)

###### WT Set #########
setwd("/Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/wt_set/pa_wt_subtractref")
temp_pa14_wt = list.files(pattern="*pa14_subtractref.txt")
all_pa14_wt_snp <- lapply(temp_pa14_wt,read.delim)
all_pa14_wt_snp
all_pa14_wt_snp_df <- as.data.frame(data.table::rbindlist(all_pa14_wt_snp,fill=T),stringsAsFactors=F,colClasses = c("character"))
all_pa14_wt_snp_df <- all_pa14_wt_snp_df[,headers] #keep headers we want
all_pa14_wt_snp_df$mutation <- paste(all_pa14_wt_snp_df$aa_ref_seq,all_pa14_wt_snp_df$aa_position,all_pa14_wt_snp_df$aa_new_seq, sep=":") #create mutation column
all_pa14_wt_snp_df$strain <- "WT_PA14"
View(all_pa14_wt_snp_df)
#add pa01 wt samples
temp_pa01_wt = list.files(pattern="*pa01_subtractref.txt")
all_pa01_wt_snp <- lapply(temp_pa01_wt,read.delim)
all_pa01_wt_snp
all_pa01_wt_snp_df <- as.data.frame(data.table::rbindlist(all_pa01_wt_snp,fill=T),stringsAsFactors=F,colClasses = c("character"))
all_pa01_wt_snp_df <- all_pa01_wt_snp_df[,headers] #keep headers we want
all_pa01_wt_snp_df$mutation <- paste(all_pa01_wt_snp_df$aa_ref_seq,all_pa01_wt_snp_df$aa_position,all_pa01_wt_snp_df$aa_new_seq, sep=":") #create mutation column
all_pa01_wt_snp_df$strain <- "WT_PA01"
View(all_pa01_wt_snp_df)
#add last batch of pa01 wt samples
setwd("/Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/wt_set/pa_wt_subtractref/wt_152")
temp_pa01_wt_152 = list.files(pattern="*_wt_pa01_subtractref.txt")
all_pa01_wt_152_snp <- lapply(temp_pa01_wt_152,read.delim)
all_pa01_wt_152_snp
all_pa01_wt_152_snp_df <- as.data.frame(data.table::rbindlist(all_pa01_wt_152_snp,fill=T),stringsAsFactors=F,colClasses = c("character"))
headers %in% colnames(all_pa01_wt_152_snp_df)  # no size column
all_pa01_wt_152_snp_df <- all_pa01_wt_152_snp_df[,c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","gene_name","gene_position","gene_product","mutation_category","new_seq","position","seq_id","snp_type","type", "title")] #keep headers we want
all_pa01_wt_152_snp_df$size <- "1"
all_pa01_wt_152_snp_df$mutation <- paste(all_pa01_wt_152_snp_df$aa_ref_seq,all_pa01_wt_152_snp_df$aa_position,all_pa01_wt_152_snp_df$aa_new_seq, sep=":") #create mutation column
all_pa01_wt_152_snp_df$strain <- "WT_PA01"
View(all_pa01_wt_152_snp_df)


all_snp_df <- rbind(all_pa01_snp_df,all_pa14_snp_df,all_pa14_wt_snp_df,all_pa01_wt_snp_df,all_pa01_wt_152_snp_df)
View(all_snp_df)
write.csv(all_snp_df,file="/Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/wt_set/all_snp_df.csv")


########## Pa01 unmapped set ########
setwd("/Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Misc_projects/Pig_wound/new_ref/pao1/unmap")
temp_pa01_unmap = list.files(pattern="*unmap_subtractref_pa01_output.txt")
pa01_unmap_snp <- lapply(temp_pa01_unmap,read.delim)
pa01_unmap_snp
pa01_unmap_snp_df <- as.data.frame(data.table::rbindlist(pa01_unmap_snp,fill=T),stringsAsFactors=F,colClasses = c("character"))
pa01_unmap_snp_df <- pa01_unmap_snp_df[,headers] #keep headers we want
pa01_unmap_snp_df$mutation <- paste(pa01_unmap_snp_df$aa_ref_seq,pa01_unmap_snp_df$aa_position,pa01_unmap_snp_df$aa_new_seq, sep=":") #create mutation column
pa01_unmap_snp_df$strain <- "unmap_Pa01"
View(pa01_unmap_snp_df)





vec_unmap <- pa01_unmap_snp_df$gene_name
length(which(table(vec_unmap)>=13)) #707
vec_unmap_13 <- names(which(table(vec_unmap)>=13))
vec_unmap_13



pa01_unmap_nocommon <- pa01_unmap_snp_df[ !(pa01_unmap_snp_df$gene_name %in% vec_unmap_13), ]
nrow(pa01_unmap_nocommon) #374
View(pa01_unmap_nocommon)


df_snp_amol_nocommon <- df_snp_amol[ !(df_snp_amol$Position %in% contig_pos_interesect), ]
nrow(df_snp_amol_nocommon) #11704
nrow(df_snp_amol) #12604
View(df_snp_amol_nocommon)
