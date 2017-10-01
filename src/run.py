import sys
import math
import traceback
from topology import *
from visualizer import *
import threading
import time

def sample():
    conformation = c.get_default_conformation()
    new_conformation = []
    for bond_length, torsion_angle in conformation:
        torsion_angle = 24.0 * (random.random()*0.5 - 0.5)
        new_conformation.append((0.5, torsion_angle))
    c.set_conformation(new_conformation)
    while True:
        new_conformation[1] = (new_conformation[1][0], new_conformation[1][1] + 0.1)
        c.set_conformation(new_conformation)
        time.sleep(0.1)

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
c = f.create_chain(villin[:3]) # just simulate the first 20
t = threading.Thread(target=sample)
t.start()

r = Visualizer(Visualizer.RES1080P, c)
r.run()
