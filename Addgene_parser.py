"""

    Plasmid parser for Addgene plasmid repository

    How does it work?
    At this stage of development you should create PlasmidParser instance with two arguments:
    id - plasmid id in Addgene plasmid repository
    path - the path on your local machine where script creates Plasmids directory.
        There will be created f'{plasmid.name}' directory with two txt files inside them:
        genebank file and csv attributes of the plasmid.

"""

from dataclasses import dataclass
from bs4 import BeautifulSoup
import requests
import pandas as pd
from urllib.request import request
import os


@dataclass()
class Description:
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


class Plasmid(Description):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def to_csv(self, path):
        if not os.path.isdir(path + f'Plasmids\\{self.name}'):
            os.makedirs(path + f'Plasmids\\{self.name}')
        with open(path + f'Plasmids\\{self.name}\\{self.name}_csv.txt', 'w', encoding='utf-8') as file:
            file.write(pd.DataFrame.from_dict({k: [v] for k, v in self.__dict__.items()}).to_csv(index_label=False))

    def to_txt(self, path, doc_seq):
        sequence = doc_seq.find_all('a', class_='genbank-file-download', href=True)[0]['href']
        req = request.Request(sequence, headers={'User-Agent': 'Mozilla/5.0'})
        seq_file = request.urlopen(req).read()

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
        # getting name
        name = self.doc.find('span', class_='material-name').text

        # getting vector backbone
        try:
            backbone = ' '.join(
                self.doc.find(string='Vector backbone').find_parent('li', class_='field').text.split()[-4::-1][-3::-1])
        except AttributeError:
            backbone = None

        # getting vector type
        try:
            vector_type = ' '.join(
                self.doc.find(string='Vector type').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            vector_type = None

        # getting selectable markers
        try:
            marker = self.doc.find(string='Selectable markers').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            marker = None

        # getting bacterial resistance(s)
        try:
            resistance = ' '.join(
                self.doc.find(string='Bacterial Resistance(s)').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            resistance = None

        # getting Growth Temperature
        try:
            growth_t = self.doc.find(string='Growth Temperature').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            growth_t = None

        # getting Growth Strain(s)
        try:
            growth_strain = self.doc.find(string='Growth Strain(s)').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            growth_strain = None

        # getting Growth instructions
        try:
            growth_instructions = ' '.join(
                self.doc.find(string='Growth instructions').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            growth_instructions = None

        # getting Copy number
        try:
            copy_num = ' '.join(self.doc.find(string='Copy number').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            copy_num = None

        # getting Gene/Insert name
        try:
            gene_insert = ' '.join(
                self.doc.find(string='Gene/Insert name').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            gene_insert = None

        # getting Total vector size (bp)
        try:
            size = int(' '.join(
                self.doc.find(string='Total vector size (bp)').find_parent('li', class_='field').text.split()[4::]))
        except AttributeError:
            size = None

        plasmid = Plasmid(name=name, gene_insert=gene_insert, growth_instructions=growth_instructions,
                          copy_num=copy_num, marker=marker,
                          growth_strain=growth_strain, resistance=resistance, vector_type=vector_type,
                          backbone=backbone, id=id,
                          vendor=self.vendor, url=self.url, growth_t=growth_t, size=size)

        plasmid.to_csv(self.path)
        plasmid.to_txt(self.path, self.doc_seq)

        PlasmidParser.plasmid_list.append(plasmid)
        return plasmid

    # only Addgene vendor is implemented yet
    def get_html(self):
        if self.vendor == 'addgene':
            url_sequence = self.url + 'sequences/'
            result_html = requests.get(self.url)
            result_seq = requests.get(url_sequence)
            return BeautifulSoup(result_html.text, 'html.parser'), BeautifulSoup(result_seq.text, 'html.parser')


def main():
    """A test function that shows how PlasmidParser works"""
    if __name__ == '__main__':
        id_list = [42888, 42876, 26248, 186478, 22222]
        PlasmidParser(id=id_list)
        plasmids_for_test = {}
        for plasmid in PlasmidParser.plasmid_list:
            plasmids_for_test.update({plasmid.name: plasmid})
        for k, v in plasmids_for_test.items():
            print(f"{k}: {v}")


main()

"""

    Plasmid parser for Addgene plasmid repository

    How does it work?
    At this stage of development you should create PlasmidParser instance with two arguments:
    id - plasmid id in Addgene plasmid repository
    path - the path on your local machine where script creates Plasmids directory.
        There will be created f'{plasmid.name}' directory with two txt files inside them:
        genebank file and csv attributes of the plasmid.

"""

from dataclasses import dataclass
from bs4 import BeautifulSoup
import requests
import pandas as pd
import urllib.request as request
import os


@dataclass()
class Description:
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


class Plasmid(Description):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def to_csv(self, path):
        if not os.path.isdir(path + f'Plasmids\\{self.name}'):
            os.makedirs(path + f'Plasmids\\{self.name}')
        with open(path + f'Plasmids\\{self.name}\\{self.name}_csv.txt', 'w', encoding='utf-8') as file:
            file.write(pd.DataFrame.from_dict({k: [v] for k, v in self.__dict__.items()}).to_csv(index_label=False))

    def to_txt(self, path, doc_seq):
        sequence = doc_seq.find_all('a', class_='genbank-file-download', href=True)[0]['href']
        req = request.Request(sequence, headers={'User-Agent': 'Mozilla/5.0'})
        seq_file = request.urlopen(req).read()

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
        # getting name
        name = self.doc.find('span', class_='material-name').text

        # getting vector backbone
        try:
            backbone = ' '.join(
                self.doc.find(string='Vector backbone').find_parent('li', class_='field').text.split()[-4::-1][-3::-1])
        except AttributeError:
            backbone = None

        # getting vector type
        try:
            vector_type = ' '.join(
                self.doc.find(string='Vector type').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            vector_type = None

        # getting selectable markers
        try:
            marker = self.doc.find(string='Selectable markers').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            marker = None

        # getting bacterial resistance(s)
        try:
            resistance = ' '.join(
                self.doc.find(string='Bacterial Resistance(s)').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            resistance = None

        # getting Growth Temperature
        try:
            growth_t = self.doc.find(string='Growth Temperature').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            growth_t = None

        # getting Growth Strain(s)
        try:
            growth_strain = self.doc.find(string='Growth Strain(s)').find_parent('li', class_='field').text.split()[2::]
        except AttributeError:
            growth_strain = None

        # getting Growth instructions
        try:
            growth_instructions = ' '.join(
                self.doc.find(string='Growth instructions').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            growth_instructions = None

        # getting Copy number
        try:
            copy_num = ' '.join(self.doc.find(string='Copy number').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            copy_num = None

        # getting Gene/Insert name
        try:
            gene_insert = ' '.join(
                self.doc.find(string='Gene/Insert name').find_parent('li', class_='field').text.split()[2::])
        except AttributeError:
            gene_insert = None

        # getting Total vector size (bp)
        try:
            size = int(' '.join(
                self.doc.find(string='Total vector size (bp)').find_parent('li', class_='field').text.split()[4::]))
        except AttributeError:
            size = None

        plasmid = Plasmid(name=name, gene_insert=gene_insert, growth_instructions=growth_instructions,
                          copy_num=copy_num, marker=marker,
                          growth_strain=growth_strain, resistance=resistance, vector_type=vector_type,
                          backbone=backbone, id=id,
                          vendor=self.vendor, url=self.url, growth_t=growth_t, size=size)

        plasmid.to_csv(self.path)
        plasmid.to_txt(self.path, self.doc_seq)

        PlasmidParser.plasmid_list.append(plasmid)
        return plasmid

    # only Addgene vendor is implemented yet
    def get_html(self):
        if self.vendor == 'addgene':
            url_sequence = self.url + 'sequences/'
            result_html = requests.get(self.url)
            result_seq = requests.get(url_sequence)
            return BeautifulSoup(result_html.text, 'html.parser'), BeautifulSoup(result_seq.text, 'html.parser')


def main():
    """A test function that shows how PlasmidParser works"""
    id_list = [42888, 42876, 26248, 186478, 22222]
    PlasmidParser(id=id_list)
    plasmids_for_test = {}
    for plasmid in PlasmidParser.plasmid_list:
        plasmids_for_test.update({plasmid.name: plasmid})
    for k, v in plasmids_for_test.items():
        print(f"{k}: {v}")


if __name__ == '__main__':
    main()
