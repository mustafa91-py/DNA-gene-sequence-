from selenium import webdriver
import time
import re

"""
a = "C:/Users/musta/Desktop/chromedriver.exe"
accession = "https://www.ncbi.nlm.nih.gov/nuccore/%s"
b = "https://www.ncbi.nlm.nih.gov/nuccore/?term=pax6"
def ncbı(title="nucletide",name="pax6"):
    b = f"https://www.ncbi.nlm.nih.gov/nuccore/?term={name}"
    return b
browser = webdriver.Chrome(a)
browser.get(ncbı(name="sox7"))
sorgu = browser.find_elements_by_css_selector("div[class=rprt]")
print(sorgu[0].text)
m = re.search("(Accession: ).*(GI)", str(sorgu[0].text))
gen_id = m.group()[len(m.group(1)):-2]
browser.get(accession % (gen_id))
time.sleep(1)
print("--"*50)
print(browser.title)
browser.close()"""


class Nucleotide(object):
    def __init__(self, name):
        self.name = name
        self.url = "https://www.ncbi.nlm.nih.gov/"
        self.browser = webdriver.Chrome("C:/Users/musta/Desktop/chromedriver.exe")
        self.browser.minimize_window()
        a = self.url + f"/nuccore/?term={self.name}"
        self.browser.get(a)
        self.cıktı = [i.text.split("\n") for i in self.browser.find_elements_by_css_selector("div[class=rprt]")]

    def accession(self):
        liste = []
        for i in self.__repr__():
            acses = re.match("(Accession: ).*(GI)", i[-3])
            k = (i[1], acses.group(0)[len(acses.group(1)):-3])
            liste.append(k)
        return liste

    def page_len(self):
        c = self.browser.find_element_by_css_selector("h3[class=page]")
        return c.text

    def __repr__(self):
        return list(self.cıktı)

    def close(self):
        self.browser.close()
    def page_next(self):
        id = 'EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Entrez_Pager.Page'
        classs = 'active page_link next'
        try:
            self.browser.find_element_by_xpath(f"//a[@id='{id}' and @class='{classs}']").click()
            self.cıktı = [i.text.split("\n") for i in self.browser.find_elements_by_css_selector("div[class=rprt]")]
        except:
            pass
    def get_fasta(self, accession):
        url = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}?report=fasta"
        self.browser.get(url)
        time.sleep(1.5)
        try:
            t = self.browser.find_element_by_tag_name("pre")
        except:
            return None
        return t.text
