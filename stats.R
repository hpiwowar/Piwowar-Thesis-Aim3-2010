# to get latex to work, first in the R GUI go to 
# Misc, start X11 server
# And execute the following line
#Sys.setenv( PATH=paste(Sys.getenv("PATH"),"/usr/texbin",sep=":") ) 

#library(Rserve)
#Rserve(args="--no-save")

#setwd("/mydir/")
source("extra_functions.R")
source("plot_summary_formula_response_CIs.R")


#### READ DATA
dat.raw = read.csv("rawdata.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(dat.raw) = gsub("_", ".", names(dat.raw))
names.pretty = cat(names(dat.raw), sep="\n")
#print names.pretty
#save(dat.raw, file="dat_raw.Rdata")
#load("dat_raw.Rdata")


#dim(dat.raw)
#summary(dat.raw)
#names(dat.raw)
#names(dat.raw)[0:20]


####  CONVERT STRINGS TO NUMBERS 
library(plyr)
dat.nums = colwise(as.numeric)(dat.raw)
library(psych)
described.dat.nums = psych::describe(dat.nums)
#write.table(described.dat.nums,"described_dat_nums.txt",append=F,quote=F,sep="\t",row.names=T)

#### START TO BUILD STATS DATASET
dat = data.frame(pmid = as.numeric(dat.nums$pmid))


# skip issn
# skip essn
# My naming convention is that "ago" means number of years between 2010 and the original variable
dat$years.ago.tr = tr(2010 - dat.nums$pubmed.year.published)
# skip date_in_pubmed
# skip authors
# Break journal into top contenders
#dat$journal.bmc.genomics                = ifelse("BMC Genomics" == dat.raw$pubmed.journal, 1, 0)
#dat$journal.cancer.res                  = ifelse("Cancer Res" == dat.raw$pubmed.journal, 1, 0)
#dat$journal.pnas                        = ifelse("Proc Natl Acad Sci U S A" == dat.raw$pubmed.journal, 1, 0)
#dat$journal.j.biol.chem                 = ifelse("J Biol Chem" == dat.raw$pubmed.journal, 1, 0)
#dat$journal.physiol.genomics            = ifelse("Physiol Genomics" == dat.raw$pubmed.journal, 1, 0)
#dat$journal.plos.one                    = ifelse("PLoS One" == dat.raw$pubmed.journal, 1, 0)

dat$num.authors.tr = tr(dat.nums$pubmed.number.authors)

####### FIRST AUTHOR VARIABLES

# skip first_author_first_name
dat$first.author.female                 = ordered(dat.nums$author.first.author.female)
dat$first.author.female[which(pmax(dat.nums$author.first.author.male, dat.nums$author.first.author.female) <= 0)] = NA

# Removing male because so highly correlated with female
#dat$first.author.male                   = ordered(dat.nums$author.first.author.male)
dat$first.author.gender.not.found       = ordered(dat.nums$author.first.author.gender.not.found)
dat$first.author.num.prev.pubs.tr       = tr(dat.nums$first.author.num.prev.pubs)
dat$first.author.num.prev.pmc.cites.tr  = tr(dat.nums$first.author.num.prev.pmc.cites)

dat$first.author.year.first.pub.ago.tr  = tr( 2010 - dat.nums$first.author.year.first.pub)
# NAs are either because there was no cluster, or the cluster included no prior papers
# (unfortunately, I didn't collect an indication of which was true, but in theory there should
# be a cluster for all papers published before 2009 Author-ity data was run, so I will 
# assume it is the latter, that there we no prior papers.)  
# I will code no prior papers as zero years since the prior paper
dat$first.author.year.first.pub.ago.tr[is.na(dat.nums$first.author.year.first.pub)] = tr(0)

dat$first.author.num.prev.microarray.creations.tr = tr(dat.nums$first.author.num.prev.microarray.creations)
dat$first.author.num.prev.oa.tr         = tr(dat.nums$first.author.num.prev.oa)
dat$first.author.num.prev.other.sharing.tr  = tr(    
                                                dat.nums$first.author.num.prev.genbank.sharing +
                                                dat.nums$first.author.num.prev.pdb.sharing +
                                                dat.nums$first.author.num.prev.swissprot.sharing)
dat$first.author.num.prev.geoae.sharing.tr = tr(dat.nums$first.author.num.prev.geoae.sharing)
dat$first.author.num.prev.geo.reuse.tr  = tr(dat.nums$first.author.num.prev.geo.reuse)
dat$first.author.num.prev.multi.center.tr = tr(dat.nums$first.author.num.prev.multi.center)
# skip first.author.prev.meta.analysis because not enough > 0

# set to NA if the number of publications is unreallistically big, because this probably means that
# Author-ity combined two or more authors together.  Better to set to NA and deal with bias than introduce
# such large numbers into our data
first.author.thresh.cluster.size = quantile(dat$first.author.num.prev.pubs.tr, .98, na.rm=TRUE)
first.author.unrealistic.pub.clusters = which(dat$first.author.num.prev.pubs.tr > first.author.thresh.cluster.size)
dat$first.author.num.prev.pubs.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.pmc.cites.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.year.first.pub.ago.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.microarray.creations.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.oa.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.other.sharing.tr [first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.geoae.sharing.tr[first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.geo.reuse.tr [first.author.unrealistic.pub.clusters] = NA
dat$first.author.num.prev.multi.center.tr[first.author.unrealistic.pub.clusters] = NA


####### LAST AUTHOR VARIABLES

# skip last_author_first_name
dat$last.author.female                 = ordered(dat.nums$last.first.author.female)
dat$last.author.female[which(pmax(dat.nums$last.first.author.male, dat.nums$last.first.author.female) <= 0)] = NA

# Removing male because so highly correlated with female
#dat$last.author.male                   = ordered(dat.nums$last.first.author.male)
dat$last.author.gender.not.found       = ordered(dat.nums$last.first.author.gender.not.found)
dat$last.author.num.prev.pubs.tr        = tr(dat.nums$last.author.num.prev.pubs)
dat$last.author.num.prev.pmc.cites.tr   = tr(dat.nums$last.author.num.prev.pmc.cites)

dat$last.author.year.first.pub.ago.tr   = tr(2010 - dat.nums$last.author.year.first.pub)
# See comment on first author above for treatment of NAs for this variable
dat$last.author.year.first.pub.ago.tr[is.na(dat.nums$last.author.year.first.pub)] = tr(0)

dat$last.author.num.prev.microarray.creations.tr = tr(dat.nums$last.author.num.prev.microarray.creations)
dat$last.author.num.prev.oa.tr          = tr(dat.nums$last.author.num.prev.oa)
dat$last.author.num.prev.other.sharing.tr   = tr(     
                                                dat.nums$last.author.num.prev.genbank.sharing +
                                                dat.nums$last.author.num.prev.pdb.sharing +
                                                dat.nums$last.author.num.prev.swissprot.sharing )
dat$last.author.num.prev.geoae.sharing.tr = tr(dat.nums$last.author.num.prev.geoae.sharing)
dat$last.author.num.prev.geo.reuse.tr   = tr(dat.nums$last.author.num.prev.geo.reuse)
dat$last.author.num.prev.multi.center.tr = tr(dat.nums$last.author.num.prev.multi.center)
# skip last.author.prev.meta.analysis because not enough > 0

# set to NA if the number of publications is unreallistically big, because this probably means that
# Author-ity combined two or more authors together.  Better to set to NA and deal with bias than introduce
# such large numbers into our data
last.author.thresh.cluster.size = quantile(dat$last.author.num.prev.pubs.tr, .98, na.rm=TRUE)
last.author.unrealistic.pub.clusters = which(dat$last.author.num.prev.pubs.tr >= last.author.thresh.cluster.size)
dat$last.author.num.prev.pubs.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.pmc.cites.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.year.first.pub.ago.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.microarray.creations.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.oa.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.other.sharing.tr [last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.geoae.sharing.tr[last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.geo.reuse.tr [last.author.unrealistic.pub.clusters] = NA
dat$last.author.num.prev.multi.center.tr[last.author.unrealistic.pub.clusters] = NA


####  SOME FUNDING VARIABLES AND STUDY VARIABLES
# skip address
# skip institution, use instition_clean instead below
# skip country, use country_clean instead below
# skip grant_numbers, use lots of other related variables below
dat$num.grant.numbers.tr = tr(dat.nums$num.grant.numbers)
dat$num.nih.is.nci.tr = tr(dat.nums$nih.is.nci)
dat$num.nih.is.nhlbi.tr = tr(dat.nums$nih.is.nhlbi)
dat$num.nih.is.ncrr.tr = tr(dat.nums$nih.is.ncrr)
dat$num.nih.is.niehs.tr = tr(dat.nums$nih.is.niehs)
dat$num.nih.is.ninds.tr = tr(dat.nums$nih.is.ninds)
dat$num.nih.is.niddk.tr = tr(dat.nums$nih.is.niddk)
dat$num.nih.is.nigms.tr = tr(dat.nums$nih.is.nigms)
dat$num.nih.is.niaid.tr = tr(dat.nums$nih.is.niaid)

# skip is_review because all 0s, as it should be
# need to use pubmed.medline.status to figure out what should really be NA because of incomplete mesh terms
dat$pubmed.is.humans    = ordered(medline.values(dat.nums$pubmed.is.humans, dat.raw$pubmed.medline.status))
dat$pubmed.is.animals   = ordered(medline.values(dat.nums$pubmed.is.animals, dat.raw$pubmed.medline.status))
dat$pubmed.is.mice      = ordered(medline.values(dat.nums$pubmed.is.mice, dat.raw$pubmed.medline.status))
dat$pubmed.is.fungi     = ordered(medline.values(dat.nums$pubmed.is.fungi, dat.raw$pubmed.medline.status))
dat$pubmed.is.bacteria  = ordered(medline.values(dat.nums$pubmed.is.bacteria, dat.raw$pubmed.medline.status))
dat$pubmed.is.plants    = ordered(medline.values(dat.nums$pubmed.is.plants, dat.raw$pubmed.medline.status))
dat$pubmed.is.viruses   = ordered(medline.values(dat.nums$pubmed.is.viruses, dat.raw$pubmed.medline.status))
dat$pubmed.is.cultured.cells    = ordered(medline.values(dat.nums$pubmed.is.cultured.cells, dat.raw$pubmed.medline.status))
# is cancer appears to depend not entirely on MEDLINE completeness
dat$pubmed.is.cancer        = ordered(dat.nums$pubmed.is.cancer)
# Open access filter, and some others below, are by journal not article, so doesn't depend on MEDLINE completion status
dat$pubmed.is.open.access   = ordered(dat.nums$pubmed.is.open.access)
dat$pubmed.is.effectiveness = ordered(medline.values(dat.nums$pubmed.is.effectiveness, dat.raw$pubmed.medline.status))
dat$pubmed.is.diagnosis     = ordered(medline.values(dat.nums$pubmed.is.diagnosis, dat.raw$pubmed.medline.status))
dat$pubmed.is.prognosis     = ordered(medline.values(dat.nums$pubmed.is.prognosis, dat.raw$pubmed.medline.status))
dat$pubmed.is.core.clinical.journal = ordered(dat.nums$pubmed.is.core.clinical.journal)
# skip is_clinical_trial because not enough
# skip is_randomized_clinical_trial because not enough
# skip is_meta_analysis because none (yay!)
dat$pubmed.is.comparative.study = ordered(medline.values(dat.nums$pubmed.is.comparative.study, dat.raw$pubmed.medline.status))
# skip multicenter because not enough (surprising?)
# skip validation because not enough (surprising?)
# skip funded_stimulus because all 0s
# external funding seems set sometimes even for 1s, so doesn't depend on MEDLINE completely
# commenting out external in favour of a combined nih below
#dat$pubmed.is.funded.nih.extramural   = ordered(dat.nums$pubmed.is.funded.nih.extramural)
dat$pubmed.is.funded.nih.intramural   = ordered(dat.nums$pubmed.is.funded.nih.intramural)
# Combining extramural, intramural, and whether the PMID had a link to funding in the NIH data sources
dat$pubmed.is.funded.nih            = ordered(pmax(dat.nums$pubmed.is.funded.nih.extramural, dat.nums$pubmed.is.funded.nih.intramural, !is.na(dat.nums$num.grants)))
dat$pubmed.is.funded.non.us.govt    = ordered(dat.nums$pubmed.is.funded.non.us.govt)
# Sharing in PDB and Swissprot very small, so will combine
dat$pubmed.is.shared.other          = ordered(pmax(dat.nums$pubmed.in.genbank, dat.nums$pubmed.in.pdb, dat.nums$pubmed.in.swissprot))
# geo.reuse is tiny but I want to use it anyway!
dat$pubmed.is.geo.reuse             = ordered(dat.nums$is.geo.reuse)
dat$pubmed.num.cites.from.pmc.tr    = tr(dat.nums$pubmed.number.times.cited.in.pmc)
dat$pubmed.num.cites.from.pmc.per.year = dat.nums$pubmed.number.times.cited.in.pmc/(2010 - dat.nums$pubmed.year.published)
# Comment these out for now, because they look wrong
#dat$found.by.highwire               = ordered(dat.nums$found.by.highwire)
#dat$found.by.scirus                 = ordered(dat.nums$found.by.highwire)
#dat$found.by.googlescholar          = ordered(dat.nums$found.by.highwire)
#dat$found.by.pmc                    = ordered(dat.nums$portal.pmids.found.by.pmc)


##### DEPENDENT VARIABLES
dat$dataset.in.geo                  = ordered(dat.nums$pubmed.in.geo)
# commenting this out in favour of a combined metric below
#dat$dataset.in.arrayexpress         = ordered(dat.nums$in.arrayexpress)
dat$dataset.in.geo.or.ae            = ordered(dat.nums$in.ae.or.geo)


###### INSTITUTION AND CORRESPONDING AUTHOR VARIABLES
dat$country.usa                    = ordered(ifelse("USA" == dat.raw$country.clean, 1, 0))
dat$country.japan                  = ordered(ifelse("Japan" == dat.raw$country.clean, 1, 0))
dat$country.germany                = ordered(ifelse("Germany" == dat.raw$country.clean, 1, 0))
dat$country.canada                 = ordered(ifelse("Canada" == dat.raw$country.clean, 1, 0))
dat$country.uk                     = ordered(ifelse("UK" == dat.raw$country.clean, 1, 0))
dat$country.china                  = ordered(ifelse("China" == dat.raw$country.clean, 1, 0))
dat$country.france                 = ordered(ifelse("France" == dat.raw$country.clean, 1, 0))
dat$country.korea                  = ordered(ifelse("Korea" == dat.raw$country.clean, 1, 0))
dat$country.australia              = ordered(ifelse("Australia" == dat.raw$country.clean, 1, 0))
dat$institution.harvard            = ordered(ifelse("Harvard University" == dat.raw$institution.clean, 1, 0))
dat$institution.nci                = ordered(ifelse("National Cancer Institute" == dat.raw$institution.clean, 1, 0))
dat$institution.stanford           = ordered(ifelse("Stanford University" == dat.raw$institution.clean, 1, 0))
# skip some institutions because not enough
#dat$institution.tokyo.daigaku      = ordered(ifelse("Tokyo Daigaku" == dat.raw$institution.clean, 1, 0))
#dat$institution.johns.hopkins      = ordered(ifelse("Johns Hopkins University" == dat.raw$institution.clean, 1, 0))
dat$institution.is.medical         = ordered(dat.nums$institution.hospital.or.medlcal)
dat$institution.rank               = dat.nums$institution.rank
dat$institution.is.higher.ed       = ordered(ifelse("Higher educ." == dat.raw$institution.sector, 1, 0))
dat$institution.is.higher.ed[which(is.na(dat$institution.rank))] = NA
dat$institution.is.govnt           = ordered(ifelse("Government" == dat.raw$institution.sector, 1, 0))
dat$institution.is.govnt[which(is.na(dat$institution.rank))] = NA
dat$institution.research.output.tr = tr(dat.nums$institution.output)
dat$institution.international.collaboration = dat.nums$international.collaboration
dat$institution.mean.norm.impact.factor     = dat.nums$institution.norm.sjr
dat$institution.mean.norm.citation.score    = dat.nums$institution.norm.citation.score

### JOURNAL VARIABLES
# skip total cites 2008 because it could be due to so many things, not very informative
#dat$journal.cites.2008.tr       = tr(dat.nums$journal.2008.cites)
dat$journal.impact.factor.log   = log.tr(dat.nums$journal.impact.factor)

# PLoS One is a big part of our sample, and doesn't have an official impact factor yet.
# Important, because it is relatively low, so better to give it our estimate instead of leave missing
# Estimate of 3 comes from http://pbeltrao.blogspot.com/2009/04/guestimating-plos-one-impact-factor.html
dat$journal.impact.factor.log[which(dat.raw$pubmed.journal== "PLoS One")] = log.tr(3)

dat$journal.5yr.impact.factor.log = log.tr(dat.nums$journal.5yr.impact.factor)
dat$journal.immediacy.index.log = log.tr(dat.nums$journal.immediacy.index)
dat$journal.num.articles.2008.tr = tr(dat.nums$journal.num.articles.2008)
dat$journal.cited.halflife      = dat.nums$journal.cited.halflife


# skip journal_policy_try because it is redundant
dat$journal.policy.requires.microarray.accession = ordered(ifelse(dat.nums$journal.policy.requires.microarray.accession > 0, 1, 0))
# says must deposit should be true whenever requires microarray accession is tru
dat$journal.policy.says.must.deposit    = ordered(ifelse(pmax(dat.nums$journal.policy.says.must.deposit, dat.nums$journal.policy.requires.microarray.accession) > 0, 1, 0))
dat$journal.policy.at.least.requests.sharing.microarray = ordered(ifelse(dat.nums$journal.policy.at.least.requests.sharing.microarray > 0, 1, 0))
dat$journal.policy.mentions.any.sharing = ordered(ifelse(dat.nums$journal.policy.mentions.any.sharing > 0, 1, 0))
dat$journal.policy.general.statement    = ordered(ifelse(dat.nums$journal.policy.general.statement > 0, 1, 0))
# requests accession should be true whenenever requires accession is true, but looks like I didn't code it that way
dat$journal.policy.requests.accession   = ordered(ifelse(pmax(dat.nums$journal.policy.requests.accession, dat.nums$journal.policy.requires.microarray.accession) > 0, 1, 0))
dat$journal.policy.mentions.exceptions  = ordered(ifelse(dat.nums$journal.policy.mentions.exceptions > 0, 1, 0))
dat$journal.policy.mentions.consequences    = ordered(ifelse(dat.nums$journal.policy.mentions.consequences > 0, 1, 0))
dat$journal.policy.contains.word.microarray = ordered(ifelse(dat.nums$journal.policy.contains.word.microarray > 0, 1, 0))
dat$journal.policy.contains.word.miame.mged = ordered(ifelse(dat.nums$journal.policy.contains.word.MIAME.MGED > 0, 1, 0))
dat$journal.policy.contains.word.arrayexpress   = ordered(ifelse(dat.nums$journal.policy.contains.word.arrayexpress > 0, 1, 0))
dat$journal.policy.contains.word.geo.omnibus    = ordered(ifelse(dat.nums$journal.policy.contains.word.geo.omnibus > 0, 1, 0))
dat$journal.policy.requests.sharing.other.data  = ordered(ifelse(dat.nums$journal.policy.requests.sharing.other.data > 0, 1, 0))
# For some reason this variable didn't have NAs, but it should have the same NAs as the other journal policy vars
dat$journal.policy.requests.sharing.other.data[is.na(dat$journal.policy.requests.accession)] = NA

# Add a variable saying how many microarray creating (as defined by being in this set)
# papers the journal has published
table.journal.counts = table(dat.raw$pubmed.journal)
dat$journal.microarray.creating.count.tr = tr(as.numeric(table.journal.counts[dat.raw$pubmed.journal]))

#### MORE FUNDING VARIABLES
# set the missing values to 0 or tr(0), as appropriate

# skip pmid:1 as duplicate column
dat$num.grants.via.nih.tr   = tr(dat.nums$num.grants)
dat$num.grants.via.nih.tr[is.na(dat$num.grants.via.nih.tr)] = tr(0)
dat$max.grant.duration.tr   = tr(dat.nums$longest.num.years)
dat$max.grant.duration.tr[is.na(dat$max.grant.duration.tr)] = tr(0)

dat$nih.first.year.ago.tr   = tr(2010 - dat.nums$first.year)
dat$nih.first.year.ago.tr[is.na(dat$nih.first.year.ago.tr)] = NA
dat$nih.last.year.ago.tr    = tr(2010 - dat.nums$last.year)
dat$nih.last.year.ago.tr[is.na(dat$nih.last.year.ago.tr)] = NA
dat$nih.cumulative.years.tr = tr(dat.nums$cumulative.years)
dat$nih.cumulative.years.tr[is.na(dat$nih.cumulative.years.tr)] = tr(0)

dat$nih.sum.avg.dollars.tr  = tr(dat.nums$sum.avg.dollars)
dat$nih.sum.avg.dollars.tr[is.na(dat$nih.sum.avg.dollars.tr)] = tr(0)
dat$nih.sum.sum.dollars.tr  = tr(dat.nums$sum.sum.dollars)
dat$nih.sum.sum.dollars.tr[is.na(dat$nih.sum.sum.dollars.tr)] = tr(0)
dat$nih.max.max.dollars.tr  = tr(dat.nums$max.max.dollars)
dat$nih.max.max.dollars.tr[is.na(dat$nih.max.max.dollars.tr)] = tr(0)


# skip
# skip grant_activity_codes list, because captured below
dat$has.R01.funding = ordered(dat.nums$has.R01.funding)
dat$has.R01.funding[is.na(dat$has.R01.funding)] = 0
#dat$has.T32.funding = ordered(dat.nums$has.T32.funding)
#dat$has.T32.funding[is.na(dat$has.T32.funding)] = 0
# other funding types dropped because don't occurr often enough

# look for patterns here
dat$has.P.funding = ordered(ifelse(grepl("P", dat.raw$grant.activity.codes), 1, 0))
dat$has.R.funding = ordered(ifelse(grepl("R", dat.raw$grant.activity.codes), 1, 0))
dat$has.T.funding = ordered(ifelse(grepl("T", dat.raw$grant.activity.codes), 1, 0))
dat$has.U.funding = ordered(ifelse(grepl("U", dat.raw$grant.activity.codes), 1, 0))
dat$has.K.funding = ordered(ifelse(grepl("K", dat.raw$grant.activity.codes), 1, 0))

# see below this block for handling of their NAs
# only pick one funding level, otherwise too singular
# not enough post2007 
dat$num.post2003.morethan500k.tr = tr(dat.nums$num.post2003.morethan500k)
dat$num.post2004.morethan500k.tr = tr(dat.nums$num.post2004.morethan500k)
dat$num.post2005.morethan500k.tr = tr(dat.nums$num.post2005.morethan500k)
dat$num.post2006.morethan500k.tr = tr(dat.nums$num.post2006.morethan500k)
#dat$num.post2007.morethan500k.tr = tr(dat.nums$num.post2007.morethan500k)
dat$num.post2003.morethan750k.tr = tr(dat.nums$num.post2003.morethan750k)
dat$num.post2004.morethan750k.tr = tr(dat.nums$num.post2004.morethan750k)
dat$num.post2005.morethan750k.tr = tr(dat.nums$num.post2005.morethan750k)
dat$num.post2006.morethan750k.tr = tr(dat.nums$num.post2006.morethan750k)
#dat$num.post2007.morethan750k.tr = tr(dat.nums$num.post2007.morethan750k)
dat$num.post2003.morethan1000k.tr = tr(dat.nums$num.post2003.morethan1000k)
dat$num.post2004.morethan1000k.tr = tr(dat.nums$num.post2004.morethan1000k)
dat$num.post2005.morethan1000k.tr = tr(dat.nums$num.post2005.morethan1000k)
dat$num.post2006.morethan1000k.tr = tr(dat.nums$num.post2006.morethan1000k)
#dat$num.post2007.morethan1000k.tr = tr(dat.nums$num.post2007.morethan1000k)

dat$num.post2003.morethan500k.tr[is.na(dat$num.post2003.morethan500k.tr)] = tr(0)
dat$num.post2004.morethan500k.tr[is.na(dat$num.post2004.morethan500k.tr)] = tr(0)
dat$num.post2005.morethan500k.tr[is.na(dat$num.post2005.morethan500k.tr)] = tr(0)
dat$num.post2006.morethan500k.tr[is.na(dat$num.post2006.morethan500k.tr)] = tr(0)
#dat$num.post2007.morethan500k.tr[is.na(dat$num.post2007.morethan500k.tr)] = tr(0)
dat$num.post2003.morethan750k.tr[is.na(dat$num.post2003.morethan750k.tr)] = tr(0)
dat$num.post2004.morethan750k.tr[is.na(dat$num.post2004.morethan750k.tr)] = tr(0)
dat$num.post2005.morethan750k.tr[is.na(dat$num.post2005.morethan750k.tr)] = tr(0)
dat$num.post2006.morethan750k.tr[is.na(dat$num.post2006.morethan750k.tr)] = tr(0)
#dat$num.post2007.morethan750k.tr[is.na(dat$num.post2007.morethan750k.tr)] = tr(0)
dat$num.post2003.morethan1000k.tr[is.na(dat$num.post2003.morethan1000k.tr)] = tr(0)
dat$num.post2004.morethan1000k.tr[is.na(dat$num.post2004.morethan1000k.tr)] = tr(0)
dat$num.post2005.morethan1000k.tr[is.na(dat$num.post2005.morethan1000k.tr)] = tr(0)
dat$num.post2006.morethan1000k.tr[is.na(dat$num.post2006.morethan1000k.tr)] = tr(0)
#dat$num.post2007.morethan1000k.tr[is.na(dat$num.post2007.morethan1000k.tr)] = tr(0)

dat$dataset.in.geo.or.ae.int = dat.nums$in.ae.or.geo

#dim(dat)

#save(dat, file="dat.Rdata")


#To do
# look into why the found by highwire values look wrong?
# remove sql for K and has U funding... the sql has a bug
# Some author variables have "author" prefix and some don't in raw... fix!
# also the name "last_first_author"... fix in original code and then here in stats file
# fix name institution.hospital.or.medlcal
# if "requires accession" then requests accession should also be true

############ Look at the data

dat.untransformed = dat

for (col in names(dat.untransformed)) {
    if (substring(col, nchar(col)-2) == ".tr") {
        print(col)
        if (substring(col, nchar(col)-6) == "log.tr") {        
            dat.untransformed[,col] = undo.log.tr(dat[,col])
        } else {
            dat.untransformed[,col] = undo.tr(dat[,col])
        }
    }
}



library(Hmisc)
library(psych)
# both Hmisc and psych have a describe
Hmisc::describe(dat.untransformed)
#psych::describe(dat)


s = summary(dataset.in.geo.or.ae.int ~ ., dat.untransformed)
s
plot.summary.formula.response.CIs(s, cex.labels=0.1, cex=0.7)
title=("Proportion of studies with datasets found in GEO or ArrayExpress")

label(dat.untransformed$first.author.num.prev.microarray.creations.tr) <- "first.author.num.prev.micro.produce.tr"
label(dat.untransformed$last.author.num.prev.microarray.creations.tr) <- "last.author.num.prev.micro.produce.tr"
label(dat.untransformed$journal.policy.at.least.requests.sharing.microarray) <- "journal.policy.requests.or.requires.micro.data"
label(dat.untransformed$journal.policy.contains.word.geo.omnibus) <- "journal.policy.contains.word.geo.omnibus"

# to see formatting, try ?dotchart2 or ?plot.summary.formula.response
mynames = names(dat.untransformed)
sections = 5
section.length = ceil(length(mynames) / sections)
#for (ii in 0:0) {
for (ii in 0:(sections-1)) {
    ## Check to make sure this is an even divide do don't miss any!
    first.count = ii*section.length
    last.count = min(length(mynames), ((ii+1)*section.length)-1)
    mynames.section = mynames[first.count:last.count]
    if ("dataset.in.geo.or.ae.int" %nin% mynames.section) {
        mynames.section = c(mynames.section, "dataset.in.geo.or.ae.int")
    }
    print(ii)
    print(mynames.section)
    s = summary(dataset.in.geo.or.ae.int ~ ., dat.untransformed[,names(dat.untransformed) %in% mynames.section], continuous = 3)
    #s
    #quartz()
    filename = paste("dotplot-vars-", ii, ".tiff", sep="")  #or pdf
    print(filename)
    tiff(filename, width=7, height=10, units="in", compression="none", res=300)
    par(mai=c(1,6,0.2,1)) # bottom, left, top, right
    plot.summary.formula.response.CIs(s, width.factor=2, cex.labels=0.4, cex=0.9, xlim=c(0,1), xlab="Proportion of studies with datasets\nfound in GEO or ArrayExpress", main="")
    dev.off()
}



#margin.table(table(dat$first.author.female, dat$dataset.in.geo.or.ae.int), 1)
#prop.table(table(dat$first.author.female, dat$dataset.in.geo.or.ae.int), 1)
#prop.table(table(dat$last.author.female, dat$dataset.in.geo.or.ae.int), 1)

######## FIGURE
library(sciplot)
par(bg="white")
lineplot.CI(x.factor = pubmed.year.published, 
            response = in.ae.or.geo, 
            data = dat.nums, 
            subset = dat.nums$pubmed.year.published<2010, 
            xlab="Year article published", ylab="Proportion of articles with datasets found in GEO or ArrayExpress")
title(main="Proportion of articles with shared datasets, by year")

tiff("figure1.tiff", width=7, height=7, units="in", compression="none", res=300)
par(bg="white")
lineplot.CI(x.factor = pubmed.year.published, 
            response = in.ae.or.geo, 
            data = dat.nums, 
            subset = dat.nums$pubmed.year.published<2010, 
            xlab="Year article published", ylab="Proportion of articles with datasets found in GEO or ArrayExpress")
title(main="Proportion of articles with shared datasets, by year")
dev.off()

############## Prep data for stats

dat.indep = dat[,!names(dat) %in% c("dataset.in.geo", "dataset.in.geo.or.ae", "dataset.in.geo.or.ae.int")]
dat.indep.stats = dat.indep[, !names(dat.indep) %in% c("pmid")]

############## Get correlation matrix

# Need to do special sorts of correlations because have binary data, not just this:
# mycor = rcorr(as.matrix(dat.indep.stats))$r

library(polycor)
myhetcorr = hetcor.modified(dat.indep.stats, use="pairwise.complete.obs", std.err=FALSE, pd=FALSE)
mycor.unadjusted = myhetcorr$correlations
#write.table(mycor.unadjusted,"/Users/hpiwowar/stats link/mycor.unadjusted.txt",append=F,quote=F,sep="\t",row.names=T)


# Some of my correlations are NAs.  Give them reasonable values
print(mycor.unadjusted["first.author.female", "first.author.gender.not.found"])
print(mycor.unadjusted["last.author.female", "last.author.gender.not.found"])

mycor.unadjusted["first.author.female", "first.author.gender.not.found"] = 0
mycor.unadjusted["first.author.gender.not.found", "first.author.female"] = 0
mycor.unadjusted["last.author.female", "last.author.gender.not.found"] = 0
mycor.unadjusted["last.author.gender.not.found", "last.author.female"] = 0

mycor.unadjusted[which(is.na(mycor.unadjusted))] = 1

# Now fix the correlation matrix if it is not positive-definite

mycor = adjust.to.positive.definite(mycor.unadjusted)

#display
library(gplots)
#quartz(height=10, width=10)
#heatmap.2(mycor, col=bluered(16), cexRow=0.5, cexCol = .8, symm = TRUE, dend = "row", trace = "none", main = "Thesis Data", margins=c(15,15), key=FALSE, keysize=0.1)

# Now, for interest, display it with data sharing variables also
#mycor.all = hetcor.modified(dat, use="pairwise.complete.obs", std.err=FALSE, pd=FALSE)
#mycor.data.sharing = mycor.all$correlations["dataset.in.geo.or.ae",]
#mycor.data.sharing.relevant = mycor.data.sharing[!names(mycor.data.sharing) %in% c("pmid", "dataset.in.geo", "dataset.in.geo.or.ae", "dataset.in.geo.or.ae.int")]
#data.sharing.colours = colorpanel(20,low="red",high="green")[10 * (1 + round(mycor.data.sharing.relevant, 1))]
#heatmap.3(mycor, ColSideColors=data.sharing.colours, col=cm.colors, cexRow=0.5, cexCol = .8, symm = TRUE, dend = "row", trace = "none", main = "Thesis Data", margins=c(15,15), key=FALSE, keysize=0.1)

#tiff("heatmap.tiff", height=10, width=10, units="in", compression="none", res=300)
#heatmap.2(mycor, col=bluered(16), cexRow=0.5, cexCol = .8, symm = TRUE, dend = "row", trace = "none", main = "Thesis Data", margins=c(15,15), key=FALSE, keysize=0.1)
#dev.off

############## FIRST ORDER ANALYSIS

##############  Determine number of First-Order factors

# Determine Number of Factors to Extract
#library(nFactors)
#eigenvectors.1st <- eigen(mycor) # get eigenvalues
# this line takes a long time
#aparallel.1st <- parallel(subject=nrow(dat.indep.stats), var=ncol(dat.indep.stats), rep=100, cent=.05)
#scree.results.1st <- nScree(eigenvectors.1st$values, aparallel.1st$eigen$qevpea)
#summary(scree.results.1st)
#plotnScree(scree.results.1st) 

# Pull out the "Optimal Coordinate"
# defined in nScree help page as 
# The optimal coordinates (OC) corresponds to an extrapolation of the preceding eigenvalue by a regression line between the eignvalue coordinates and the last eigenvalue coordinate
#number.factors.1st = scree.results.1st$Components$noc
number.factors.1st = 15


##############  Do First-Order Factor Analysis

# Maximum liklihood doesn't converge because too 
#fit.ml = factanal(dat, number.factors.1st, rotation="promax", covmat=mycor)

# Use princip axis when maximum liklihood fails to converge:
library(psych)
fit.fa.1st = fa(mycor, number.factors.1st, fm="minres", rotate="promax", 
                scores=FALSE, residuals=TRUE, n.obs=max(dim(dat.indep.stats)))

#to show the loadings sorted by absolute value
print(fit.fa.1st, sort=TRUE)

#fa.diagram(fit.fa.1st)
#cluster.plot(fit.fa.1st)

# look at residuals
#fit.fa.1st.residuals = fit.fa.1st$residual
#heatmap.2(fit.fa.1st.residuals, col=bluered(16), cexRow=0.5, cexCol = .8, symm = TRUE, dend = "row", trace = "none", main = "Thesis Data", margins=c(15,15), key=FALSE, keysize=0.1)


factor.names.1st = c(
"MR1"="Large NIH grant",
"MR2"="Has journal policy",
"MR3"="NOT institution NCI or intramural",
"MR4"="Count of R01 & other NIH grants",
"MR5"="Journal impact",
"MR6"="Last author num prev pubs & first year pub",
"MR7"="Journal policy consequences & long halflife",
"MR8"="Institution high citations & collaboration",
"MR9"="NOT geo reuse & YES high institution output",
"MR10"="NOT animals or mice",
"MR11"="Humans & cancer",
"MR12"="Instititution is government & NOT higher ed",
"MR13"="NOT K funding or P funding",
"MR14"="Authors prev GEOAE sharing & OA & microarray creation",
"MR15"="First author num prev pubs & first year pub")

for (afactor in names(factor.names.1st)) {
    #print(paste(afactor, ": ", factor.names.1st[afactor], sep=""))
    print(paste(factor.names.1st[afactor], sep=""))
    print.thresh(fit.fa.1st$loadings[, afactor], .4, TRUE)
    print.thresh(fit.fa.1st$loadings[, afactor], -0.4, FALSE)
}



############## SECOND ORDER ANALYSIS

##############  Determine number of Second-Order factors
# Equations as per Factor Analysis, Second Edition, by Gorsuch

# Correlation of my first order results
fit.fa.1st.cor = fit.fa.1st$r.scores

#eigenvectors.2nd <- eigen(fit.fa.1st.cor) # get eigenvalues
#aparallel.2nd <- parallel(subject=nrow(fit.fa.1st.cor), var=ncol(fit.fa.1st.cor), rep=100, cent=.05)
#scree.results.2nd <- nScree(eigenvectors.2nd$values, aparallel.2nd$eigen$qevpea)
#scree.results.2nd
#plotnScree(scree.results.2nd) 

#number.factors.2nd = scree.results.2nd$Components$noc
number.factors.2nd = 5


##############  Do Second-Order Factor Analysis

# Ideally uncorrelated, but want it to be a good fit
#fit.fa.2nd = fa(fit.fa.1st.cor, number.factors.2nd, fa="minres", rotate="varimax")

fit.fa.2nd = fa(fit.fa.1st.cor, number.factors.2nd, fa="minres", rotate="varimax")
print(fit.fa.2nd, sort=TRUE)

#fa.diagram(fit.fa.2nd)
#cluster.plot(fit.fa.2nd)

##############  Map variables directly to second order analysis
# Equations as per Factor Analysis, Second Edition, by Gorsuch

U = fit.fa.2nd$uniquenesses * diag(number.factors.1st)
P = fit.fa.2nd$loadings
A = cbind(P, U)
Pvf = fit.fa.1st$loadings

Pvo = Pvf %*% A

############# HAVE A LOOK AT THESE RESULTS

#Pvo[1:24, 1:4]
# Interesting:  last.author.num.prev.geoae.sharing.tr          0.134139157 -0.066308705  0.138307260 -0.04026715
# what about how it would correlate with our main endpoint?


factor.names.2nd = c(
"MR1"="Amount of NIH funding",
"MR2"="Cancer & humans",
"MR3"="OA journal & previous GEO-AE sharing",
"MR4"="Journal impact factor and policy",
"MR5"="Higher Ed in USA")

# On original variables
for (afactor in names(factor.names.2nd)) {
    #print(paste(afactor, ": ", factor.names.2nd[afactor], sep=""))    
    print(paste(factor.names.2nd[afactor], sep=""))    
    print.thresh(Pvo[, afactor], .3, TRUE)
    print.thresh(Pvo[, afactor], -0.3, FALSE)
}

rownames(fit.fa.2nd$loadings) <- factor.names.1st
# On first-order factors
for (afactor in names(factor.names.2nd)) {
#    print(paste(afactor, ": ", factor.names.2nd[afactor], sep=""))
    print(paste(factor.names.2nd[afactor], sep=""))
    print.thresh(fit.fa.2nd$loadings[, afactor], .3, TRUE)
    print.thresh(fit.fa.2nd$loadings[, afactor], -0.3, FALSE)
}
########   IMPUTE VARIABLES SO CAN CALCULATE SCORES


# Just want one output variables
dat$dataset.in.geo.or.ae.int = dat.nums$in.ae.or.geo
dat.impute.input = dat[,!names(dat) %in% c("dataset.in.geo", "dataset.in.geo.or.ae")]


# Show the pattern of NAs
library(mice)
#md.pattern(dat.impute.input)

### Impute variables
#mice.output = mice(dat.impute.input)  # singular
mice.output = mice(dat.impute.input[,c(2:109, 114:117, 122)], m=1, maxit=3)

 
# Now flush out the rest of the scores 
dat.imputed = complete(mice.output, 1)
dat.for.scores = dat.imputed[,!names(dat.imputed) %in% c("dataset.in.geo.or.ae.int")]
dat.for.scores$num.post2003.morethan500k.tr = dat.indep.stats$num.post2003.morethan500k.tr
dat.for.scores$num.post2004.morethan500k.tr = dat.indep.stats$num.post2004.morethan500k.tr
dat.for.scores$num.post2005.morethan500k.tr = dat.indep.stats$num.post2005.morethan500k.tr
dat.for.scores$num.post2006.morethan500k.tr = dat.indep.stats$num.post2006.morethan500k.tr
dat.for.scores$num.post2003.morethan750k.tr = dat.indep.stats$num.post2003.morethan750k.tr
dat.for.scores$num.post2004.morethan750k.tr = dat.indep.stats$num.post2004.morethan750k.tr
dat.for.scores$num.post2005.morethan750k.tr = dat.indep.stats$num.post2005.morethan750k.tr
dat.for.scores$num.post2006.morethan750k.tr = dat.indep.stats$num.post2006.morethan750k.tr
dat.for.scores$num.post2003.morethan1000k.tr = dat.indep.stats$num.post2003.morethan1000k.tr
dat.for.scores$num.post2004.morethan1000k.tr = dat.indep.stats$num.post2004.morethan1000k.tr
dat.for.scores$num.post2005.morethan1000k.tr = dat.indep.stats$num.post2005.morethan1000k.tr
dat.for.scores$num.post2006.morethan1000k.tr = dat.indep.stats$num.post2006.morethan1000k.tr

###### COMPUTE FIRST ORDER SCORES

scores.1st = factor.scores.bartlett(dat.for.scores, fit.fa.1st)


####### FIRST ORDER REGRESSION
# Not sure if this will be a primary result

dat.regress = data.frame(dataset.in.geo.or.ae.int = dat$dataset.in.geo.or.ae.int)
dat.regress[,dimnames(scores.1st)[[2]]] = scores.1st

# If I name them here, then setting up the regressions is more awkward due to long names
#names(dat.regress) = c("dataset.in.geo.or.ae.int", factor.names.1st[names(dat.regress)[-1]])

library(rms)
dd.regress = datadist(dat.regress)
options(datadist='dd.regress')
options(digits=2)

# Remove MR13, MR1, MR3
# Remove interactions for MR3, MR8, MR12
f.1st.nonlinear.interactions.reduced = lrm(formula = dataset.in.geo.or.ae.int ~ 
    (
    rcs(MR2, 3) +
        rcs(MR4, 3) + rcs(MR14, 3) + 
        rcs(MR13, 3) +
        rcs(MR5, 3) + rcs(MR7, 3)  + rcs(MR8, 3) + rcs(MR10, 3) + rcs(MR12, 3) + rcs(MR6, 3) + 
        rcs(MR1, 3) + rcs(MR1, 3) + rcs(MR11, 3) + 
        rcs(MR9, 3) + rcs(MR15, 3) )^2,
    dat=dat.regress, x=T, y=T)
anova(f.1st.nonlinear.interactions.reduced)

summ.1st.nonlinear.interactions.reduced = summary(f.1st.nonlinear.interactions.reduced)
summ.1st.nonlinear.interactions.reduced.dimnames = dimnames(summ.1st.nonlinear.interactions.reduced)[[1]][seq(1,length(dimnames(summ.1st.nonlinear.interactions.reduced)[[1]]),2)]
dimnames(summ.1st.nonlinear.interactions.reduced)[[1]][seq(1,length(dimnames(summ.1st.nonlinear.interactions.reduced)[[1]]),2)] = factor.names.1st[summ.1st.nonlinear.interactions.reduced.dimnames]
par(bg="white")
plot.summary.rms.norangelabels(summ.1st.nonlinear.interactions.reduced, q = c(0.95), col=gray(0.5), log=T, cex=0.9, width=0.01, cex.c=0.001, at=c(0.25, 0.5, 1, 2, 4, 8))
title("Multivariate nonlinear regressions with interactions")

### Dots of first-order factors
dat.regress.named = dat.regress
names(dat.regress.named) = c("dataset.in.geo.or.ae.int", factor.names.1st[names(dat.regress)[-1]])

dots.1st.nonlinear.interactions.reduced = summary(dataset.in.geo.or.ae.int ~ .,
    dat=dat.regress.named)
plot.summary.formula.response.CIs(dots.1st.nonlinear.interactions.reduced, width.factor=2, cex.labels=0.5, cex=0.7, 
    xlab="Percentage of studies with datasets found in GEO or ArrayExpress", 
    main="Univariate data sharing behaviour on first order factors")
#plot.summary.formula.response
#?summary.formula

tiff("dotplot-firstorder.tiff", width=7, height=7, units="in", compression="none", res=300)
par(mai=c(1,6,0.1,1)) # bottom, left, top, right
plot.summary.formula.response.CIs(dots.1st.nonlinear.interactions.reduced, width.factor=2, cex.labels=0.4, cex=0.9, xlim=c(0,1), 
    xlab="Proportion of studies with datasets\nfound in GEO or ArrayExpress", 
    main="")
#    main="Univariate data sharing behaviour on first order factors")
dev.off()

			
######### SECOND ORDER REGRESSION


loadings.2nd = Pvo

#fit.pa.2nd.tovars = list(loadings=loadings.2nd[,(1+length(colnames(fit.pa.2nd$weights))):length(colnames(loadings.2nd))])
fit.pa.2nd.tovars = list(loadings=loadings.2nd[,colnames(fit.fa.2nd$loadings)])
fit.pa.2nd.tovars$uniquenesses = apply(fit.pa.2nd.tovars$loadings^2, 1, sum)

scores.to.dat.2nd = factor.scores.bartlett(dat.for.scores, fit.pa.2nd.tovars)

dat.regress.2nd = data.frame(dataset.in.geo.or.ae.int = dat$dataset.in.geo.or.ae.int) 
dat.regress.2nd[,dimnames(scores.to.dat.2nd)[[2]]] = scores.to.dat.2nd

library(rms)

dd.regress.2nd = datadist(dat.regress.2nd)
options(datadist='dd.regress.2nd')
options(digits=2)

f.2nd = lrm(formula = dataset.in.geo.or.ae.int ~ ., dat=dat.regress.2nd, x=T, y=T)
anova(f.2nd)

### Now try with nonlinear splines:
f.2nd.nonlinear = lrm(formula = dataset.in.geo.or.ae.int ~ 
    rcs(MR1, 4) + rcs(MR2, 4) + rcs(MR3, 4) + rcs(MR4, 4) + rcs(MR5, 4),
    dat=dat.regress.2nd, x=T, y=T)
anova(f.2nd.nonlinear)

f.2nd.nonlinear.interactions = lrm(formula = dataset.in.geo.or.ae.int ~ 
    (rcs(MR1, 4) + rcs(MR2, 4) + rcs(MR3, 4) + rcs(MR4, 4) + rcs(MR5, 4))^2,
    dat=dat.regress.2nd, x=T, y=T)
anova(f.2nd.nonlinear.interactions)

f.2nd.nonlinear.interactions.reduced = lrm(formula = dataset.in.geo.or.ae.int ~ 
    (rcs(MR3, 3) + rcs(MR1, 3) + rcs(MR4, 3) + rcs(MR5, 3)  + rcs(MR2, 3))^2,
    dat=dat.regress.2nd, x=T, y=T)
anova(f.2nd.nonlinear.interactions.reduced)

summ.2nd.nonlinear.interactions.reduced = with(dat.regress.2nd, summary(f.2nd.nonlinear.interactions.reduced))
summ.2nd.nonlinear.interactions.reduced
summ.2nd.nonlinear.interactions.reduced.dimnames = dimnames(summ.2nd.nonlinear.interactions.reduced)[[1]][seq(1,length(dimnames(summ.2nd.nonlinear.interactions.reduced)[[1]]),2)]
dimnames(summ.2nd.nonlinear.interactions.reduced)[[1]][seq(1,length(dimnames(summ.2nd.nonlinear.interactions.reduced)[[1]]),2)] = factor.names.2nd[summ.2nd.nonlinear.interactions.reduced.dimnames]
par(bg="white")
plot.summary.rms.norangelabels(summ.2nd.nonlinear.interactions.reduced, q = c(0.95), col=gray(0.5), log=T, cex=0.9, width=0.01, cex.c=0.001, at=c(0.25, 0.5, 1, 2, 4))
title("Multivariate nonlinear regression with interactions")

plot.summary.rms.norangelabels(summ.2nd.nonlinear.interactions.reduced, q = c(0.95), col=gray(0.5), log=T, cex=1, width=0.01)

### Dots of second-order factors
dat.regress.2nd.named = dat.regress.2nd
names(dat.regress.2nd.named) = c("dataset.in.geo.or.ae.int", factor.names.2nd[names(dat.regress.2nd)[-1]])

dots.2nd.nonlinear.interactions.reduced = summary(dataset.in.geo.or.ae.int ~ .,
    dat=dat.regress.2nd.named)
dots.2nd.nonlinear.interactions.reduced
plot.summary.formula.response.CIs(dots.2nd.nonlinear.interactions.reduced, width.factor=2, cex.labels=0.5, cex=0.7)

tiff("dotplot-secondorder.tiff", width=7, height=5, units="in", compression="none", res=300)
par(mai=c(1,6,0.1,1)) # bottom, left, top, right
plot.summary.formula.response.CIs(dots.2nd.nonlinear.interactions.reduced, width.factor=2, cex.labels=0.4, cex=0.9, xlim=c(0,1), 
    xlab="Proportion of studies with datasets\nfound in GEO or ArrayExpress", 
    main="")
#    main="Univariate data sharing behaviour on second order factors")
dev.off()


#########  CLUSTERING

if (FALSE) {
    n <- nrow(dat.regress.2nd[,2:6])
    wss <- rep(0, 20)
    wss[1] <- (n - 1) * sum(apply(dat.regress.2nd[,2:6], 2, var))
    for (i in 2:20) {
    	wss[i] <- sum(kmeans(dat.regress.2nd[,2:6], centers = i)$withinss)
    }
    plot(1:20, wss, type = "b", xlab = "Number of groups", ylab = "Within groups sum of squares")

    k = 8
    km = kmeans(dat.regress.2nd[,2:6], centers=k)
    dat.explore.clusters.2nd = cbind(dat, dat.regress.2nd, km.cluster.2nd = km$cluster)
    cbind(aggregate(dat.explore.clusters.2nd[,125:130], by=list(kmeans(dat.regress.2nd[,2:6], centers=k)$cluster), mean), table(km$cluster))


    library(clValid)
    dat.cluster.2nd = dat.regress.2nd[,2:6]
    dat.cluster.2nd.tmp = dat.regress.2nd[1:100,2:6]
    stab = clValid(dat.cluster.2nd.tmp, 2:3, clMethods = c("hierarchical", "kmeans", "pam"), validation = "stability")
    optimalScores(stab)

    op <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    plot(stab, measure = c("APN", "AD", "ADM"), legend = FALSE)
    plot(nClusters(stab), measures(stab, "APN")[, , 1], type = "n", axes = F, xlab = "", ylab = "")
    legend("center", clusterMethods(stab), col = 1:9, lty = 1:9, pch = paste(1:9))
    par(op)

    stab.internal = clValid(dat.cluster.2nd.tmp, 2:3, clMethods = c("hierarchical", "kmeans", "pam"), validation = "internal")
    summary(stab.internal)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    plot(stab.internal, legend = FALSE)
    plot(nClusters(stab.internal), measures(stab.internal, "Dunn")[, , 1], type = "n", axes = F, xlab = "", ylab = "")
    legend("center", clusterMethods(stab.internal), col = 1:9, lty = 1:9, pch = paste(1:9))
    par(op)
}

library(tree)
g.smallest = tree(dataset.in.geo.or.ae.int ~ 
    rcs(MR1, 4) + rcs(MR2, 4) + rcs(MR3, 4) + rcs(MR4, 4) + rcs(MR5, 4),
    dat = dat.regress.2nd,
    control=tree.control(nobs=max(dim(dat.regress.2nd)), minsize=5))

g.bigger = tree(dataset.in.geo.or.ae.int ~ 
    rcs(MR1, 4) + rcs(MR2, 4) + rcs(MR3, 4) + rcs(MR4, 4) + rcs(MR5, 4),
    dat = dat.regress.2nd,
    control=tree.control(nobs=max(dim(dat.regress.2nd)), mindev=0.003))
    
g = g.smallest    
g
tree.screens()
plot(g, type="u")
text(g, digits=2)
pred.tree = predict(g)
pred.full.grouped = cut2(dat.regress.2nd$dataset.in.geo.or.ae.int, seq(0, 1, by=0.01)) 
tile.tree(g, pred.full.grouped)   
close.screen(all = TRUE)
    
# Just split on sharing and cancer, at medians:
# I think this is my take-away:
#aggregate(dat.regress.2nd$dataset.in.geo.or.ae.int, by=list(cut2(dat.regress.2nd$MR3, g=2), cut2(dat.regress.2nd$MR2, g=2)), mean)

#counts
numerator = aggregate(dat.regress.2nd$dataset.in.geo.or.ae.int, by=list(cut2(dat.regress.2nd$MR3, g=2), cut2(dat.regress.2nd$MR2, g=2)), function(x) sum(x==1))
denominator = aggregate(dat.regress.2nd$dataset.in.geo.or.ae.int, by=list(cut2(dat.regress.2nd$MR3, g=2), cut2(dat.regress.2nd$MR2, g=2)), function(x) sum(x<2))

binconf(numerator$x, denominator$x, method="wilson") 


aa = c(1,2)
bb = c(3,4)
cc = c(1,3)
dd = c(2,4)
combos = cbind(aa, bb, cc, dd)
for (ii in 1:4) {
    indices = combos[,ii]
    print(indices)
    print(numerator$x[indices])
    print(sum(numerator$x[indices]))
    print(sum(denominator$x[indices]))
    print(binconf(sum(numerator$x[indices]), sum(denominator$x[indices]), method="wilson") )
}

indices = 1:4
print(numerator$x[indices])
print(sum(numerator$x[indices]))
print(sum(denominator$x[indices]))
print(binconf(sum(numerator$x[indices]), sum(denominator$x[indices]), method="wilson") )


#write.table(dat.explore.clusters.2nd,"/Users/hpiwowar/stats link/at_explore_clusters_2nd.txt",append=F,quote=F,sep="\t",row.names=T)

#dim(dat.explore.clusters.2nd)
#names(dat.explore.clusters.2nd)



########

#save.image("image.RData")
#setwd("/Users/hpiwowar/Documents/Code/hpiwowar/pypub/trunk/src/aim3/stats")
#load("image.RData")

#setwd("/Users/hpiwowar/Documents/Code/hpiwowar/pypub/trunk/src/aim3/stats")
#for (iii in 1:4) source("aim3_stats_20100217b.R", echo=TRUE)

dat.journals = dat.raw
dat.journals$dataset.in.geo.or.ae.int = as.numeric(dat.journals$in.ae.or.geo)

# Consolidate the less frequent journals into Other
top.journals = names(sort(table(dat.journals$pubmed.journal), decreasing=TRUE)[0:50])
dat.journals$journal.factor = factor(dat.journals$pubmed.journal)
levels(dat.journals$journal.factor)[levels(dat.journals$journal.factor) %nin% top.journals] <- "OTHER"


### This is a good plot
s = summary(dataset.in.geo.or.ae.int ~ journal.factor, dat.journals)
s
plot.summary.formula.response.CIs(s, cex.labels=0.5, xlim=c(0,1))

filename = paste("dotplot-journals.tiff", sep="")
print(filename)
tiff(filename, width=7, height=10, compression="none", units="inches", res=300)  # or can do pdf
par(mai=c(1,6,0.2,1)) # bottom, left, top, right
plot.summary.formula.response.CIs(s, width.factor=2, cex.labels=0.5, cex=0.9, xlim=c(0,1), xlab="Proportion of studies with datasets\nfound in GEO or ArrayExpress", main="")
dev.off()



dat.untransformed$journal.factor = dat.journals$journal.factor

s = summary(dataset.in.geo.or.ae.int ~ journal.cited.halflife, dat.journals, dat.journals)
plot.summary.formula.response.CIs(s, cex.labels=0.5, xlim=c(0,1))

boxplot(journal.cited.halflife ~ dataset.in.geo.or.ae.int, dat=dat.untransformed)
table(dat.journals[which(as.numeric(dat.journals$journal.cited.halflife)<2),]$pubmed.journal)

####################

# GRAPHS for presentation


quartz()
par(omi=c(1.5, 0, 0, 0))
par(bg="white")
bargraph.CI.ordered(x.factor = journal.factor, 
            response = dataset.in.geo.or.ae.int, 
            data = dat.journals, 
            err.width=0.03,
            col=grey(0.7), 
            ylim=c(0,1),
            cex=0.8,
            cex.names=0.7,
            las=2, 
            space=c(0, .5),
            ylab="Proportion of datasets shared", 
            err.col=grey(.5))

dat.institutions = dat.raw
dat.institutions$dataset.in.geo.or.ae.int = as.numeric(dat.institutions$in.ae.or.geo)
dat.institutions$institutions.factor = factor(dat.institutions$institution.clean)
top.institutions = names(sort(table(dat.institutions$institutions.factor), decreasing=TRUE)[0:25])
display.institutions = c(top.institutions, "University of Pittsburgh", "The University of British Columbia")
levels(dat.institutions$institutions.factor)[levels(dat.institutions$institutions.factor) %nin% display.institutions] <- "OTHER"

quartz()
par(omi=c(2, 0, 0, 0))
par(bg="white")
bargraph.CI.ordered(x.factor = institutions.factor, 
            response = dataset.in.geo.or.ae.int, 
            data = dat.institutions, 
            err.width=0.03,
            col=grey(0.7), 
            ylim=c(0,1),
            cex=0.8,
            cex.names=0.9,
            las=2, 
            space=c(0, .5),
            ylab="Proportion of datasets shared", 
            err.col=grey(.5))

dat.funding = dat.raw
dat.sum.sum.dollars = cut(dat.funding)
quartz()
par(omi=c(0, 0, 0, 0))
par(bg="white")
bargraph.CI(x.factor = cut(dat.institutions$institution.rank, breaks=seq(0, 2000, 100)), 
            response = dataset.in.geo.or.ae.int, 
            data = dat.institutions, 
            err.width=0.03,
            col=grey(0.7), 
            ylim=c(0,1),
            cex=0.8,
            cex.names=0.9,
            las=2, 
            space=c(0, .5),
            names.arg = seq(1, 2000, 100),
            ylab="Proportion of datasets shared", 
            err.col=grey(.5))
#bargraph.CI(x.factor = cut(log(dat.institutions$journal.impact.factor), breaks=seq(0, 4, 0.5)), 
            

library(gplots)
quartz(height=5, width=5)
row.names(fit.fa.1st.cor) = factor.names.1st[dimnames(fit.fa.1st$loadings)[[2]]]
dimnames(fit.fa.1st.cor)[[2]] = factor.names.1st[dimnames(fit.fa.1st$loadings)[[2]]]
par(bg="white")
heatmap.2(fit.fa.1st.cor, col=bluered(16), cexRow=0.8, cexCol=0.8, symm = TRUE, 
    dend = "row", trace = "none", main = "Thesis Data", margins=c(15,15), key=FALSE, keysize=0.1)

