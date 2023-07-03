"""

    Plasmid parser for Addgene plasmid repository

    How does it work?
    At this stage of development you should create PlasmidParser instance with two arguments:
    id - plasmid id in Addgene plasmid repository
    path - the path on your local machine where script creates Plasmids directory.
        There will be created f'{plasmid.name}' directory with two txt files inside them:
        genebank file and csv attributes of the plasmid.

"""
import json
from dataclasses import dataclass
from bs4 import BeautifulSoup
import requests
import pandas as pd
import urllib.request as request
import os


@dataclass()
class Description:
    """A class that contains all the information about the plasmid that is given on Addgene site.
       This class represents a plasmid and provides attributes to store various information about the plasmid, 
       such as its name, identifier, vendor, URL, size, backbone, vector type, marker, resistance, growth temperature, 
       growth strain, growth instructions, copy number, and gene insertion name.

       :param name: The name of the plasmid that is given on Addgene site
       :type name: str
       :param id: Vendor's plasmid identifier
       :type id: int
       :param vendor: Which vendor's site to parse
       :type vendor: str
       :param url: URL to the plasmid (e.g. https://www.addgene.org/22222/)
       :type url: str
       :param size: Size of the plasmid in base pairs (bp)
       :type size: int
       :param backbone: The ancestor, vector backbone of the plasmid
       :type backbone: str
       :param vector_type: The purpose of the plasmid (e.g. mammalian expression, bacterial expression)
       :type vector_type: list of str
       :param marker: The type of selectable markers plasmid has (e.g. Gentamicin)
       :type marker: str
       :param resistance: Types of bacterial resistance (e.g. Kanamycin, 50 μg/mL)
       :type resistance: str
       :param growth_t: Growth temperature to grow in bacteria (e.g. 37°C)
       :type growth_t: str
       :param growth_strain: Which strain to use to grow this plasmid in bacteria (e.g. DH5alpha, xl1-blue)
       :type growth_strain: str
       :param growth_instructions: Some specified information about growth in bacteria
       :type growth_instructions: str
       :param copy_num: The characteristic of the plasmid: the number of copies of a given plasmid in a cell
       :type copy_num: str
       :param gene_insert: The insertion/gene name according to the authors
       :type gene_insert: str
       :param sequence: Transformed gbk file into text; a sequence with annotations
       :type sequence: str"""
    name: str  # The name of the plasmid that is given on Addgene site
    id: int  # Vendor's plasmid identificator
    vendor: str  # Which vendor's site to parse
    url: str  # URL to the plasmid (e.g. https://www.addgene.org/22222/)
    size: int  # Size of the plasmid in base pairs (bp)
    backbone: str  # The ancestor, vector backbone of the plasmid
    vector_type: list  # The purpose of the plasmid (e.g. mammalian expression, bacterial expression)
    marker: str  # The type of selectable markers plasmid has (e.g. Gentamicin)
    resistance: str  # Types of bacterial resistance (e.g. Kanamycin, 50 μg/mL)
    growth_t: str  # Growth temperature to grow in bacteria (e.g. 37°C)
    growth_strain: str  # Which strain use to grow this plasmid in bacteria (e.g. DH5alpha, xl1-blue)
    growth_instructions: str  # Some specified information about growth in bacteria
    copy_num: str  # The characteristic of the plasmid: the number of copies of a given plasmid in a cell.
    gene_insert: str  # The insertion/gene name according to the authors
    sequence: str  # Transformed gbk file into text; a sequence with annotations


class Plasmid(Description):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def to_csv(self, path):
        if not os.path.isdir(path + f'Plasmids\\{self.name}'):
            os.makedirs(path + f'Plasmids\\{self.name}')
        with open(path + f'Plasmids\\{self.name}\\{self.name}_csv.txt', 'w', encoding='utf-8') as file:
            file.write(pd.DataFrame.from_dict({k: [v] for k, v in self.__dict__.items()}).to_csv(index_label=False))

    def to_json(self, path):
        if not os.path.isdir(path + f'Plasmids\\{self.name}'):
            os.makedirs(path + f'Plasmids\\{self.name}')
        with open(path + f'Plasmids\\{self.name}\\{self.name}.json', 'w', encoding='utf-8') as file:
            file.write(json.dumps({k: [v] for k, v in self.__dict__.items()}))

    def to_txt(self, path, seq_file):
        if not os.path.isdir(path + f'Plasmids\\{self.name}'):
            os.makedirs(path + f'Plasmids\\{self.name}')
        with open(path + f'Plasmids\\{self.name}\\{self.name}.txt', 'wb') as file:
            file.write(seq_file)

        # if there is no info about total vector size on the parsed page
        if self.size is None:
            with open(path + f'Plasmids\\{self.name}\\{self.name}.txt', 'r') as file:
                self.size = int(file.readline().split()[2])


class PlasmidParser:
    plasmid_list = []

    def __init__(self, id, base_url: str = "https://www.addgene.org/", id_start: int = id, id_end: int = None,
                 vendor: str = 'addgene',
                 path: str = f''):
        self.id = id
        self.base_url = base_url
        self.id_start = id_start
        self.id_end = id_end
        self.vendor = vendor
        self.path = path
        if isinstance(self.id, int):
            self.url = self.base_url + f'{self.id}/'
            self.doc, self.doc_seq = self.get_html()

            self.get(self.id)

        if isinstance(self.id, list):
            for plasmid_id in self.id:
                self.url = self.base_url + f'{plasmid_id}/'
                self.doc, self.doc_seq = self.get_html()

                self.get(plasmid_id)

    def get(self, id: int):
        if self.sorry_defence():
            return None
        else:
            plasmid = Plasmid(name=self.get_name(), gene_insert=self.get_insert(),
                              growth_instructions=self.get_growth_instructions(),
                              copy_num=self.get_copy_number(), marker=self.get_marker(),
                              growth_strain=self.get_growth_strain(), resistance=self.get_resistance(),
                              vector_type=self.get_resistance(),
                              backbone=self.get_backbone(), id=id,
                              vendor=self.vendor, url=self.url, growth_t=self.get_growth_t(), size=self.get_size(),
                              sequence=self.get_txt(self.doc_seq))

            # plasmid.to_csv(self.path) # Uncomment if you want to write down a text with csv
            # plasmid.to_json(self.path) # Uncomment if you want to write down a json file
            plasmid.to_txt(self.path, self.get_txt(self.doc_seq))

            PlasmidParser.plasmid_list.append(plasmid)
            return plasmid

    # only Addgene vendor is implemented yet
    def get_html(self):
        if self.vendor == 'addgene':
            url_sequence = self.url + 'sequences/'
            result_html = requests.get(self.url)
            result_seq = requests.get(url_sequence)
            return BeautifulSoup(result_html.text, 'html.parser'), BeautifulSoup(result_seq.text, 'html.parser')

    def get_txt(self, doc_seq):
        sequence = doc_seq.find_all('a', class_='genbank-file-download', href=True)[0]['href']
        req = request.Request(sequence, headers={'User-Agent': 'Mozilla/5.0'})
        seq_file = request.urlopen(req).read()
        return seq_file

    def sorry_defence(self):
        if self.vendor == 'addgene':
            sorry = self.doc.find(string='Sorry!')
            if sorry == 'Sorry!':
                print(f"Sorry! There is no such ID.")
                return True
            else:
                return False

    def get_name(self):
        # getting name
        name = self.doc.find('span', class_='material-name').text
        return name

    def get_backbone(self):
        # getting vector backbone
        try:
            backbone = ' '.join(
                self.doc.find(string='Vector backbone').find_parent('li', class_='field').text.split()[-4::-1][-3::-1])
        except AttributeError:
            backbone = None
        return backbone

    def get_vector_type(self):
        # getting vector type
        try:
            vector_type = ' '.join(
                self.doc.find(string='Vector type').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            vector_type = None
        return vector_type

    def get_marker(self):
        # getting selectable markers
        try:
            marker = ' '.join(
                self.doc.find(string='Selectable markers').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            marker = None
        return marker

    def get_resistance(self):
        # getting bacterial resistance(s)
        try:
            resistance = ' '.join(
                self.doc.find(string='Bacterial Resistance(s)').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            resistance = None
        return resistance

    def get_growth_t(self):
        # getting Growth Temperature
        try:
            growth_t = ' '.join(
                self.doc.find(string='Growth Temperature').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            growth_t = None
        return growth_t

    def get_growth_strain(self):
        # getting Growth Strain(s)
        try:
            growth_strain = ' '.join(
                self.doc.find(string='Growth Strain(s)').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            growth_strain = None
        return growth_strain

    def get_growth_instructions(self):
        # getting Growth instructions
        try:
            growth_instructions = ' '.join(
                self.doc.find(string='Growth instructions').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            growth_instructions = None
        return growth_instructions

    def get_copy_number(self):
        # getting Copy number
        try:
            copy_num = ' '.join(self.doc.find(string='Copy number').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            copy_num = None
        return copy_num

    def get_insert(self):
        # getting Gene/Insert name
        try:
            gene_insert = ' '.join(
                self.doc.find(string='Gene/Insert name').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            gene_insert = None
        return gene_insert

    def get_size(self):
        # getting Total vector size (bp)
        try:
            size = int(' '.join(
                self.doc.find(string='Total vector size (bp)').find_parent('li', class_='field').text.split()[4::]))
        except AttributeError:
            size = None
        return size


def main():
    """A test function that shows how PlasmidParser works"""
    id_list = [1, 42888, 42876, 26248, 186478, 22222]
    PlasmidParser(id=id_list)
    plasmids_for_test = {}
    for plasmid in PlasmidParser.plasmid_list:
        plasmids_for_test.update({plasmid.name: plasmid})
    for k, v in plasmids_for_test.items():
        print(f"{k}: {v}")


if __name__ == '__main__':
    main()
