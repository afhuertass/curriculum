er #!/usr/bin/python
# -*- coding: utf-8 -

import json
import io
import pandas

import os

import operator
from nltk.corpus import stopwords

from nltk import bigrams
import string
from collections import Counter

import matplotlib.pyplot as plt

from nltk.tokenize import TweetTokenizer

tknizer = TweetTokenizer()


class TwittAnalyser():
    def __init__(self , jsonfile):
        #constructor  - twitt analyser class
        self.jsonfile = jsonfile
        self.allTokens = []
        self.allDates = []
        punctu = list(string.punctuation)
        self.stop = stopwords.words('english') + punctu + ['via' , 'rt']
        
    def lowerTokens(self, tokens):
        tokens = [ token.lower() for token in tokens]
        return tokens
    
    def AllTokens(self):
        #abrimos el archivo
        co = 0
        with io.open( self.jsonfile , "r" , encoding="utf-8-sig") as df:
            for line in df:
                ## cada linea es parseada en un objeto json
                tweet = json.loads(line)
                try:
                    text = tweet['text']
                except KeyError:
                    print("error reading key")
                    continue
                tokens = tknizer.tokenize( text  )
                tokens = self.lowerTokens( tokens )
                self.allTokens.append(tokens) # una lista de tokens
                self.allDates.append( tweet['created_at'] )
                print("tweets procs:"  + str( co ) )
                co = co +1
                #if co == 10000:
                 #   return
        tweets_total = len(self.allTokens)
        print("Tokens Sucessfully loaded, tweets: " + str(tweets_total) )
        print( str( len( self.allTokens ) ) ) 

    def FrecuencyCount(self , n = 5 ):
        counter_all = Counter()
        for tokens in self.allTokens:
            terms = [ token for token in tokens if token not in self.stop  ]
            counter_all.update(terms)

        for term in  counter_all.most_common(n):
            print term
    def getBigrams(self, n = 5):
        #
        counter_bi = Counter()
        for tokens in self.allTokens:
            terms = [ token for token in tokens if token not in self.stop  ]
            #terms = [ token for token in tokens if token not in self.stop and token.startswith( ('#','@' ) ) ]
            terms_bigrams = bigrams(terms )
            counter_bi.update(terms_bigrams)

        for bi in counter_bi.most_common(n):
            print bi

    def hashtagSeries(self,  hashtag = "#none"):
        dates = []

        for tokens in self.allTokens:
            if hashtag in tokens:
                dates.append( self.allDates[ self.allTokens.index(tokens)  ] )

        ones = [1]*len(dates)
        idx = pandas.DatetimeIndex(dates)

        HashS = pandas.Series( ones , index = idx)
        #HashS.plot.line()
        per_minute = HashS.resample('1Min').sum().fillna(0)
        per_minute.plot.line( legend = True , label= hashtag)
        
        #plt.show()


    def simplyClassifier(self ,  lexicon_pos = "" , lexicon_neg = "" , outdir =""):
        # build the simple classifier
        lex_pos = [] # lista para el vocabulario positivo
        lex_neg = []
        if not os.path.isdir( outdir ) :
            print("Error, no es directorio")
            return 
        with io.open(lexicon_pos , "r" , encoding="latin-1") as pos:
            for lex in pos:
                if lex == ";":
                    continue
                lex_pos.append( lex.split('\n')[0] )
                
        with io.open(lexicon_neg , "r" , encoding="latin-1") as neg:
            for lex in neg:
                
                if lex == ";" or lex == " " :
                    continue
                lex_neg.append( lex.split('\n')[0] )
        print(lex_neg)
        
        # listo tenemos la lista de caracteres positivos y negativos
        pos_dir = outdir+"pos/"
        neg_dir = outdir+"neg/"
        neu_dir = outdir+"neu/"
        if not os.path.exists( pos_dir):
            os.makedirs(pos_dir)
        if not os.path.exists( neg_dir ):
            os.makedirs( neg_dir )

        if not os.path.exists( neu_dir ):
            os.makedirs ( neu_dir  )
        ## ahora a clasificar un tweet por conteo lexico
        s = " "
        pos_count = 0
        neg_count = 0
        neu_count = 0
        for tokens in self.allTokens:
            score = 0
            for token in tokens:
                
                    
                if token in lex_pos:
                    score = score + 1
                if token in lex_neg:
                    
                    score = score - 1
            text = s.join( tokens )

            #print(str(score))
            if score > 0 : # clasificado como positivo
                pos_count = pos_count + 1
                pos_name =  'pos_' + str(pos_count)
                with io.open( pos_dir+pos_name , 'w' ) as posf:
                    posf.write(text)

            if score < 0 :
                neg_count = neg_count + 1
                neg_name = 'neg_' + str(neg_count)
                with io.open( neg_dir+neg_name , 'w') as negf:
                    negf.write(text)

            if score == 0:
                neu_count = neu_count + 1
                neu_name = 'neu_' + str(neu_count)
                with io.open( neu_dir+neu_name , 'w') as negf:
                    negf.write(text)
                

        
            
            
print("this is just a test")

twa = TwittAnalyser("../data/walking-dead-30-8-oct.json")
twa.AllTokens()
twa.FrecuencyCount(20)
twa.getBigrams(10)

#twa.hashtagSeries(hashtag="rick")


twa.hashtagSeries(hashtag="negan")
twa.hashtagSeries(hashtag="kill")
twa.hashtagSeries(hashtag="carol")
twa.hashtagSeries(hashtag="ezekiel")
plt.xlabel("Hora(UTC)")
plt.ylabel("Tweets por minuto ")
plt.show()


#twa.simplyClassifier( "../data/lexicon/opinion-en/positive-words.txt" , "../data/lexicon/opinion-en/negative-words.txt" , "./classi"  )

"""            
def lowerTokens(tokens):
    tokens = [ token.lower() for token in tokens]
    return tokens



def AllTokens(jsonFile = ""):
    # esta funcion retornara una lista de todos los tokens
    with io.open( tw)

twitts = "./airline-data.json"


li = 1
with  io.open( twitts , "r",encoding="utf8" ) as f:
    for line in f.readlines():
        if not line:
            print("Something went wrong, bro")
            break
        
        #print( line )
        tweet = json.loads(line)
        
        text = tweet['text']
        text2 = word_tokenize( tweet['text'])
        
        text3 =  tknzr.tokenize( tweet['text'] )
        #text2 = tweet['text'].decode('utf8')
        #print (line)
        print( li )
        print( text )
        print( text3 )
        li = li +1
            
            
       
        
        #print( text2 ) 
    
"""
