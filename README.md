protein-tools
============
Tools to generate unfolded protein models from amino-acid sequences

### Features
* Convert amino-acid sequence into unfolded topology
* Set backbone conformation using torsion angles

![Screenshot](https://github.com/saliksyed/protein-tools/blob/master/docs/screenshot.png)

### Usage

NOTE: This is an unfinished, personal, experimental project. Master may be
completely broken at any given time! This should *NOT* be used for
any academic/research work. 

```
import sys
import math
import traceback
from topology import *
from visualizer import *
import threading
import time


f = ForceField()
villin = """MTKLSAQVKGSLNITTPGLQIWRIEAMQMVPVPSSTFGSFFDGDCYIILAIHKTASSLSY
DIHYWIGQDSSLDEQGAAAIYTTQMDDFLKGRAVQHREVQGNESEAFRGYFKQGLVIRKG
GVASGMKHVETNSYDVQRLLHVKGKRNVVAGEVEMSWKSFNRGDVFLLDLGKLIIQWNGP
ESTRMERLRGMTLAKEIRDQERGGRTYVGVVDGENELASPKLMEVMNHVLGKRRELKAAV
PDTVVEPALKAALKLYHVSDSEGNLVVREVATRPLTQDLLSHEDCYILDQGGLKIYVWKG
KKANEQEKKGAMSHALNFIKAKQYPPSTQVEVQNDGAESAVFQQLFQKWTASNRTSGLGK
THTVGSVAKVEQVKFDATSMHVKPQVAAQQKMVDDGSGEVQVWRIENLELVPVDSKWLGH
FYGGDCYLLLYTYLIGEKQHYLLYVWQGSQASQDEITASAYQAVILDQKYNGEPVQIRVP
MGKEPPHLMSIFKGRMVVYQGGTSRTNNLETGPSTRLFQVQGTGANNTKAFEVPARANFL
NSNDVFVLKTQSCCYLWCGKGCSGDEREMAKMVADTISRTEKQVVVEGQEPANFWMALGG
KAPYANTKRLQEENLVITPRLFECSNKTGRFLATEIPDFNQDDLEEDDVFLLDVWDQVFF
WIGKHANEEEKKAAATTAQEYLKTHPSGRDPETPIIVVKQGHEPPTFTGWFLAWDPFKWS
NTKSYEDLKAELGNSRDWSQITAEVTSPKVDVFNANSNLSSGPLPIFPLEQLVNKPVEEL
PEGVDPSRKEEHLSIEDFTQAFGMTPAAFSALPRWKQQNLKKEKGLF"""
c = f.create_chain(villin[:40]) # just simulate the first 40 residues

# Sample Conformations:
# count = 0
# start = time.time()
# while True:
#     c.set_conformation(c.get_random_conformation())
#     # TODO: c.get_energy()
#     count += 1
#     if count%10 == 0:
#         print "%d samples evaluated in %.2f seconds" % (count, (time.time() - start))

r = Visualizer(Visualizer.RES1080P, c)
r.run()

```


