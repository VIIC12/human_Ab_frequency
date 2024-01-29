# README

Automates the search for the [Observed Antibody Space database (OAS)](https://opig.stats.ox.ac.uk/webapps/oas/). <br/>
Download the pre-filtered database via the website, it is not necessary to unpack the files. <br/>
You can then filter for specific HeavyChains (HC) according to various search criteria. <br/>

**Usage:** ./humanab.py [-h] [--IGHV IGHV] [--IGHD IGHD] [--IGHJ IGHJ] [--CDRH3_length CDRH3_LENGTH] [--h3_motif H3_MOTIF] [--database DATABASE] [--full_results FULL_RESULTS] [--outputdir OUTPUTDIR] [--overwrite OVERWRITE] <br/>

| **Option**     | **Explanation**                                                   | **Default** |
|----------------|-------------------------------------------------------------------|-------------|
| --help (-h)    | show this help message and exit                                   |             |
| --IGHV         | "3-" or "3-22" or multiple by using | like: "3-20|3-22"           |             |
| --IGHD         | "3-" or "3-22" or multiple by using | like: "3-20|3-22"           |             |
| --IGHJ         | "5" or multiple by using | like: "4|5"                            |             |
| --CDRH3_length | Length of CDHR3 region as int (WITH C-X-W, so if necessary add 2) |             |
| --h3_motif     | "." for one, ".*" for 0-many, like "YY.D.*G"                      |             |
| --database     | Set database directory                                            |             |
| --full_results | Safe full results table. 1 for True, 0 for False                  | 1           |
| --outputdir    | Output directory                                                  | "output/"   |
| --overwrite    | 1 for True, 0 for False                                           | 0           |

**Example:** ./humanab.py --IGHV "1-69" <br/> <br/>


#### <u>Results</u>
- results.csv and (if --full_results=1) fullresults.csv <br/>
- fullresults.csv is the combined database table according to your search criteria. <br/>
- results.csv looks like this: <br/>

 **null**           | **Subject1** | **...** | **Meta**                                              
--------------------|--------------|---------|-------------------------------------------------------
 Hits               | 100          |         | IGHV: 1-69, IGHD: , IGHJ: , CDRH3-length: , H3-motif: 
 Total sequences    | 1000         |         |                                                       
 Percentage of hits | 10.0         |         |                                                       
 Hits per 1 million | 100000.0     |         |                                                       


---
version = 1.0