import os

class TCLcalls:
    """A collection of functions designed to execute tcl scripts."""
    def __init__(self):
        pass

    @staticmethod
    def executeGetRotationalMatrices(psf: str, dcd: str, outputFile: str, refpdb: str) -> None:
        """"""
        with open("getRotationalMatrices.sh", 'w+') as out:
            out.write(f"#!/bin/sh\n"
                      f"alias vmd=/Applications/VMD_1.9.4a57-arm64-Rev12.app/Contents/Resources/VMD.app/Contents/MacOS/VMD\n"
                      f"vmd -dispdev text -e /Users/danielsantos/Downloads/CRALBP_analysis/scripts/getRotationalMatrices.tcl -args {outputFile} {psf} {dcd} {refpdb}")

        os.system("chmod +x getRotationalMatrices.sh")
        os.system("./getRotationalMatrices.sh")

    @staticmethod
    def executeRotateProteinToPreferred(pdb, yaw, pitch, roll, outfile) -> None:
        """"""
        with open("rotateProteinToPreferred.sh", 'w+') as out:
            out.write(f"#!/bin/sh\n"
                      f"alias vmd=/Applications/VMD_1.9.4a57-arm64-Rev12.app/Contents/Resources/VMD.app/Contents/MacOS/VMD\n"
                      f"vmd -dispdev text -e /Users/danielsantos/Downloads/CRALBP_analysis/scripts/rotateProteinToPreferred.tcl -args {pdb} {yaw} {pitch} {roll} {outfile}")

        os.system("chmod +x rotateProteinToPreferred.sh")
        os.system("./rotateProteinToPreferred.sh")


    @staticmethod
    def executeGetLipidProtInterTXT(psf, dcd, outfile) -> None:
        """"""
        with open("getLipidProtInter-txt.sh", "w+") as out:
            out.write(f"#!/bin/sh\n"
                      f"alias vmd=/Applications/VMD_1.9.4a57-arm64-Rev12.app/Contents/Resources/VMD.app/Contents/MacOS/VMD\n"
                      f"vmd -dispdev text -e /Users/danielsantos/Downloads/CRALBP_analysis/scripts/getLipidProtInter-txt.tcl -args {psf} {dcd} {outfile}")

        os.system("chmod +x getLipidProtInter-txt.sh")
        os.system("./getLipidProtInter-txt.sh")

