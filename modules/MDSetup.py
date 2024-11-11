from MDAnalysis import Universe

class MDSetup:
    """Collection of Setup functions for Molecular Dynamics Analysis"""

    def __init__(self, *, universe = None, DCD = None, PSF = None, nptboxfile = None, trajectory_length = 0):

        self._universe = universe
        self._dcd = DCD
        self._psf = PSF
        self._nptboxfile = nptboxfile
        self._trajectory_length = trajectory_length

        self.initializeUni(self._universe, self._dcd, self._psf)

    def setDCD(self, dcd):
        """Sets dcd file for current simulation object."""
        self._dcd = dcd

    def setPSF(self, psf):
        """Sets psf file for current simulation object"""
        self._psf = psf

    def setBox(self, boxfile):
        """Sets box size file for current simulation object"""
        self._nptboxfile = boxfile

    def setTrajLen(self, length):
        """Sets frame count for current simulation object"""
        self._trajectory_length = length

    def initializeUni(self, universe: None, dcd, psf):
        """Initializes the universe object."""
        if not universe:
            uni = Universe(psf, dcd)
            self._universe = uni

    def getUniverse(self):
        """Returns universe of current simulation object"""
        return self._universe

    def getDCD(self):
        """Return dcd file for current simulation object."""
        return self._dcd

    def getPSF(self):
        """Return psf file for current simulation object."""
        return self._psf

    def getTrajLen(self):
        """Returns frame count for current simulation object"""
        return self._trajectory_length



