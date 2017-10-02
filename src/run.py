import sys
import math
import traceback
from topology import *
from visualizer import *
import threading
import time

<<<<<<< HEAD
=======
def sample():
    conformation = c.get_default_conformation()
    new_conformation = []
    for bond_length, torsion_angle in conformation:
        torsion_angle = 24.0 * (random.random()*0.5 - 0.5)
        new_conformation.append((8.0, torsion_angle))
    c.set_conformation(new_conformation)
    while True:
        new_conformation[0] = (new_conformation[0][0], new_conformation[0][1] + 0.1)
        c.set_conformation(new_conformation)
        time.sleep(0.1)
>>>>>>> ae320ce2323b9c5816fe723633e5adcb09aadb63

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
