import time
class Timer(object):
    def __init__(self):
        self.elapsed={}
        self.start={}
    def tic(self, tag):
        self.start[tag] = time.time()
    def toc(self, tag):
        tend = time.time()
        assert tag in self.start, f"{tag} must be inside self.start!"
        tElapsed = tend - self.start[tag]
        if tag in self.elapsed:
          self.elapsed[tag] += tElapsed
        else:
          self.elapsed[tag] = tElapsed
    def get_elapsed_time(self,tag):
        return self.elapsed[tag]
    def dump(self):
        print("===== Timing =====")
        for tag in self.elapsed.keys():
          print(f"{tag} took {self.elapsed[tag]} s")
