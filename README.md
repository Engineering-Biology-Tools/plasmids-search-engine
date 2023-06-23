# plasmids-search-engine
Browser for plasmids search.

**Purpose**: this repositiry is the codebase for plasmids search engine.  

**What is plasmid**: according to [Wikepedia](https://en.wikipedia.org/wiki/Plasmid), plasmid is a small, extrachromosomal DNA molecule within a cell that is physically separated from chromosomal DNA and can replicate independently.

**Data sources**: browser engine store data from [Addgene](https://www.addgene.org/), [Registry of Standard Biological Parts](http://parts.igem.org/Main_Page), and other sources.  

**Technological Stack**: Python, PostgreSQL, Django, Flask, FastAPI, ...   


# Linux (Nix or NixOs)

To run the script:
```
nix-shell
python 'Addgene_parser.py'
```
