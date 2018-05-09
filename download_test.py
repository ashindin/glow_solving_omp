import requests
import json
import os

def get_direct_url(url_yandex):
    r = requests.get("https://cloud-api.yandex.net/v1/disk/public/resources/download?public_key="+url_yandex)
    if r.status_code==200:
        return json.loads(r.content)['href']
    else:
        return 1

urls_filename="./download_test.url"
fid=open(urls_filename,'r')
urls=fid.readlines()
fid.close()

for i in range(len(urls)):
    if urls[i][-1]=='\n':
        urls[i]=urls[i][:-1]

direct_urls=[]
for i in range(len(urls)):
    direct_url=get_direct_url(urls[i])
    if direct_url==1:
        print("Can not get direct url for " + urls[i])
        break
    direct_urls.append(direct_url)
#print(direct_urls)

if os.path.exists("./test.fits")==True:
    print("File test.fits already exists")
else:
    print("Downloading file test.fits")
    os.system('wget -O test.fits' +' "'+ direct_urls[0]+'"')


