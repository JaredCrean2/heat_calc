#!/usr/bin/env python3
import requests

# downloads data from NRELs PSM3 data set as a csv

api_key="3MDrHEBBQuRawvH7VPXCziYNsMNFwtFS1O9CUy2P"
url = "http://developer.nrel.gov/api/nsrdb/v2/solar/psm3-download.json?api_key={}".format(api_key)

print("url = ", url)


#payload = "names=2012&interval=60&email=jcrean01@gmail.com&wkt=POINT(35.0844, 106.6504)"

#payload="names=2012&interval=60&email=jcrean01%40gmail.com&wkt=POINT%2835.0844%2C 106.6504%29"

#payload = "names=2012&leap_day=false&interval=60&utc=false&email=jcrean01%40gmail.com&affiliation=NREL&attributes=dhi%2Cdni%2Cwind_speed%2Cair_temperature&wkt=wkt=POINT(179.9901 -16.96)"

payload = { 'names': '2012', 'interval': '60', 'email': 'jcrean01@gmail.com', 'wkt': 'Point(-106.6204, 35.0844)'}

headers = {
            'content-type': "application/x-www-form-urlencoded",
            'cache-control': "no-cache"
}

response = requests.get(url, params=payload, headers=headers)

print(response.text)

