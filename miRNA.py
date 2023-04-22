from sequence import Sequence
class MiRNA(Sequence):
    def __init__(self, h, s):
        Sequence.__init__(self, h, s)


    def getSeedSequence(self, begin, end):
        '''
        extract the seed sequence. parameters are one-based positions
        '''
        return self._sequence[begin-1:end]