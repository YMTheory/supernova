from abc import abstractmethod, ABC

class reaction_XS(ABC):

    def __init__(self, name):
        self.name = name
   
    @abstractmethod
    def diffXS(self, Ev, T, flavour):
        pass

    @abstractmethod
    def totXS(self, Ev, flavour):
        pass
