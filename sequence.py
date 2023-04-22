class Sequence:

    def __init__(self, h, s):
        '''
        Constructor
        '''
        self._header = h
        self._sequence = s
        self._gcPercent = 0.0
        self._ntCount  = 0
        self._gcCount = 0


    def getGCPercent(self):
        return self._gcPercent

    def getHeaderLine(self):
        return self._header

    def getSequence(self):
        return self._sequence


    def calcGC(self):
        '''
            calculate the gcPercent for this sequence
        '''
        self._gcCount = 0.0;
        self._ntCount = 0.0;
        n = 0
        while n < len(self._sequence):
            nt = self._sequence[n]

            if nt == 'G' or nt == 'C':
                self._gcCount += 1
                self._ntCount += 1

            if nt == 'A' or nt == 'T':
                self._ntCount +=1

            n += 1

        self._gcPercent = self._gcCount/self._ntCount







