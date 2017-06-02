#!/usr/bin/python



#miner of the data.
import time
import tweepy

from tweepy import OAuthHandler
from tweepy import Stream
from tweepy.streaming import StreamListener
   

class MyListener( StreamListener ):

    def __init__(self, filename ):
        self.outfile = filename
        
    
    def on_data(self, data):
        #do something with the data
        try:
            with open(self.outfile , 'a' ) as f:
                print("dataaa")
                f.write(str(data))
                return True
        except BaseException as e:
            print("Error on_data (waiting one minute): %s" , str(e))
            time.sleep(1*60) # dormir
        return True
            
    def on_error(self, status_code):
        #handle error
        if( status_code == 40):
            time.sleep(10*60) # parar 10 minutos
        return True
    

#for status in tweepy.Cursor(api.home_timeline).items(10):

    #print( status.text )

consumer_key = 	"QCpxppGqyu6OFecJQ3NAyy6cB"
consumer_secret = "1zZzHbXDnHsQrkDNn4nxjp78f704L4bEaotTT8IwJZiaau0ZbP"

access_token = "750504233921830912-QOhExGPYmELYN077NJY8IZ49vIZtGP3"
access_secret = "Df6v1kgb2li2eFFOiUPWeqyenXPlEjHqDTSd5aAMo4Dxw"

auth = OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_secret)
 
api = tweepy.API(auth)

    
# programa main 
companies = ['@Avianca' , '@LATAM_CO', '@VivaColombiaco' , '@LufthansaLatina' ]
print (companies)
twitter_stream = Stream( auth , MyListener("./airline-data.json" ) )
twitter_stream.filter( track = companies )

