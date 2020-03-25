
rm(list=ls())
setwd("C:/Users/mat460/Dropbox/3rd project/Real data")
#load("dataall.RDATA") 
library(tm)
library(wordcloud)
rmv.nonASCII <- function(dat){
  # convert string to vector of words
  dat2 <- unlist(strsplit(as.character(dat[1]), split=" "))
  # find indices of words with non-ASCII characters
  dat3 <- grep("dat2", iconv(dat2, "latin1", "ASCII", sub="dat2"))
  # subset original vector of words to exclude words with non-ASCII char
  if (length(dat3)!=0) {
    dat4 <- dat2[-dat3];
    # convert vector back to a string
    dat <- as.character(paste(dat4, collapse = " "));
  }
  return(dat)
}

stemCompletion_mod <- function(x,dict) {
  PlainTextDocument(stripWhitespace(paste(stemCompletion(unlist(strsplit(as.character(x)," ")),dictionary=dict, type="shortest"),sep="", collapse=" ")))
}

#dat <- read.table("iran_news2.txt", sep="\t")

text2dm <- function(dat, rmwords=c("iran","iranian","tehran"), ind=c(1:267)){
  
  dat <- sapply(dat,function(x) rmv.nonASCII(x))
  dat <- sapply(dat , function(x) gsub("(f|ht)tp(s?)://(.*)[.][a-z]+", "", x))
  # make corpus
  mycorpus <- Corpus(VectorSource(dat[ind])) # from 16 June until first of August
  #writeLines(as.character(mycorpus[[2]]))
  #writeLines(as.character(mycorpus[[2]])) 
  #stripWhitespace: eliminate extra white-spaces
  mycorpus = tm_map(mycorpus, stripWhitespace) 
  #tolower: convert text to lower case
  mycorpus = tm_map(mycorpus, content_transformer(tolower))
  #removeWords: remove words like stopwords
  mycorpus = tm_map(mycorpus, removeWords, stopwords("english"))
  #removePunctuation: remove punctuation symbols
  mycorpus = tm_map(mycorpus, removePunctuation)
  #removes the words specified in rmwords
  dict <- mycorpus <- tm_map(mycorpus, removeWords, rmwords)  
  mycorpus <- tm_map(mycorpus, stemDocument)
  #mycorpus <- tm_map(mycorpus, stemCompletion_mod, dictionary=mycorpus) 
  #mycorpus <- tm_map(mycorpus, content_transformer(function(x, dict)
  #  paste(stemCompletion(strsplit(stemDocument(x), ' ')[[1]], dict), collapse = ' ')), dict)
  tdm = TermDocumentMatrix(mycorpus)
  m = as.matrix(tdm) 
  word_freqs = sort(rowSums(m), decreasing=TRUE) 
  dm = data.frame(word=names(word_freqs), freq=word_freqs)
  return(list(dm=dm,m=m))
} 

#data= read.table("StatusID, Date,UserID,Username,text(Desired Time).txt", sep="\t")
#save(file="alldata.RDATA",data)
#load("alldata.RDATA")
ind.uni.text <- duplicated(data$V4)

data1 <- data[!ind.uni.text,]
data1 <- data1[order(data1$V3),]
dates <- unique(data1$V3)
a=numeric()
inds=list()

    
for (i in 1:length(dates)){
inds[[i]]=which(data1$V3==dates[i])
a=append(a,length(inds[[i]]))
}

data.ts1=data.ts2=data.ts3=data.ts4=character()
# t0=day one: data1$V3[1]= Ð 2009-05-05 Ð 15 ordibehesht 88
# t1= Ð 2009-10-22 Ð    30 Mehr 88
# t2=Ð 2010-03-16 Ð  25 esfand 88
# t3=last day= Ð 2010-08-08 Ð   17 mordad 89
t1=which(data1$V3==dates[145])[1]
t2=which(data1$V3==dates[290])[1]                                                                                                                                                                                                                                                                                  =which(data1$V3==dates[290])[1]
t3=which(data1$V3==dates[435])[1]
#t4=which(data1$V3==dates[400])[1]

for (i in 1:2000){  
#data.ts1[i]=paste(data1$V4[(1+(i-1)*floor(t1/2000)):(i*floor(t1/2000))],collapse = " ")
data.ts1[i]=paste(data1$V4[sample((1+(i-1)*floor(t1/2000)):(i*floor(t1/2000)), size=50, replace = FALSE)],collapse = " ")
data.ts2[i]=paste(data1$V4[(t1+1+(i-1)*floor((t2-t1)/2000)):(t1+i*floor((t2-t1)/2000))],collapse = " ")
data.ts3[i]=paste(data1$V4[(t2+1+(i-1)*floor((t3-t2)/2000)):(t2+i*floor((t3-t2)/2000))],collapse = " ")
#data.ts4[i]=paste(data1$V4[(t3+1+(i-1)*floor(t4/2000)):(t3+i*floor(t4/2000))],collapse = " ")
} 
t0=Sys.time()
rm(list=ls()[!(ls() %in% c('data.ts1',"data.ts1","data.ts2","data.t3"))])
out1 <- text2dm(data.ts1, rmwords=c("iran","iranian","tehran"), ind=c(1:2000))
out2 <- text2dm(data.ts2, rmwords=c("iran","iranian","tehran"), ind=c(1:2000)) 
out3 <- text2dm(data.ts3, rmwords=c("iran","iranian","tehran"), ind=c(1:2000)) 
#save(file="allstagetext.RDATA", out1,out2,out3,data.ts1,data.ts2,data.ts3)


#wordcloud(as.numeric(sort(a3, decreasing = TRUE)[1:20]), names(sort(a3, decreasing = TRUE)[1:20]), random.order=FALSE, colors=brewer.pal(8, "Dark2"))
#png("MachineLearningCloud.png", width=12, height=8, units="in", res=300)
#wordcloud(dm$word[5:50], dm$freq[5:50], random.order=FALSE, colors=brewer.pal(8, "Dark2"))









##################################################################################
# TWITTER ############################################################################
##################################################################################
#RESTART R session!

library(devtools)
install_github("twitteR", username="geoffjentry")
library(twitteR)

api_key <- "KBPBY4VUlYSXq2YH1aVrOyS5A"

api_secret <- 	"ByIO4jQkHL6QJxMhRPeucr8nlXWNFpMTSPZvd6DqdoHIWGqDxJ"

access_token <- 	"1674335910-mDe4PBjUc9RtPo72miQZODsmOFpeZYzGLmPqtg3"

access_token_secret <-	"Fb09myhkDIqdsbMmDPHTMpUTkej6u4IrCo3lbPwbWvIxV"

setup_twitter_oauth(api_key,api_secret,access_token,access_token_secret)

tweets <- userTimeline("NIHforHealth", n=3000)
n=length(tweets)
# collect tweets in english containing 'data mining'
start.date <- '2016-01-01'
end.date <- as.character(as.Date(start.date) + 100)
tweets <- searchTwitter("#blacklivesmatter", n=1000, lang="en",since=start.date, until=end.date)
ebola <- twListToDF(ebolaTweets)
tweets = searchTwitter("#blacklivesmatter", lang="en", n=1000)
length(tweets)

# convert tweets into a data frame
tweets_df = twListToDF(tweets)

write.table(tweets_df, "mydata.txt", sep="\t")

# extract the text content of the first tweet
tweets[[1]]$getText()

# extract the user name of the first tweet
tweets[[1]]$getScreenName()

# extract the tweet Id of the first tweet
tweets[[1]]$getId()

# extract date and time of publication of the first tweet
tweets[[1]]$getCreated()

# extract source user agent
tweets[[1]]$getStatusSource()



# extract the text content of the all the tweets
sapply(tweets, function(x) x$getText())

# extract the user name of the all the tweets
sapply(tweets, function(x) x$getScreenName())

# extract the Id number of the all the tweets
sapply(tweets, function(x) x$getId())

# extract the date and time of publication of the all the tweets
sapply(tweets, function(x) x$getCreated())

# extract the source user agent of the all the tweets
sapply(tweets, function(x) x$getStatusSource())




