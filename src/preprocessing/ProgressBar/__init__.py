__author__ = 'gwy15'

import math
import time
from multiprocessing.dummy import Lock as ThreadLock

class ProgressBar:
    global charset
    charset = "▏▎▍▌▋▊▉█"

    def __init__(self, maxCount, maxLength=50, printCount=True, printPercentage=True, printTime=False):
        self.lock = ThreadLock()
        self.maxCount = maxCount
        self.currentCount = 0
        self.maxLength = maxLength
        self.printCount = printCount
        self.printPercentage = printPercentage
        self.printTime = printTime
        self.startTime = time.time()

    def begin(self):
        self.update(0)

    def update(self, currentCount, maxCount=None):
        self.currentCount = currentCount
        if maxCount:
            self.maxCount = maxCount
        self.lock.acquire()
        self._print()
        self.lock.release()
        return
    
    def grow(self, growCount=1):
        return self.update(self.currentCount+growCount)

    def _print(self):
        percentage = float(self.currentCount)/float(self.maxCount)
        percentage = max(0., percentage)
        percentage = min(1., percentage)
        length = float(self.maxLength)*percentage
        intPart = math.floor(length)
        restPart = length-intPart

        string =  charset[-1]*intPart # print main part
        if restPart > 0.01:
            string += charset[math.floor(restPart*len(charset))] # print rest part
        string += ' ' * (self.maxLength - len(string)) # fill with blank
        string += charset[0] # print boarder

        if self.printPercentage:
            head = "%3.1f"%(percentage*100.,)
            head = (5 - len(head))*' ' + head
            string = head + "% " + string
        if self.printCount:
            countStr = str(self.currentCount)+'/'+str(self.maxCount)
            string += countStr
        if self.printTime:
            secondsSpent = math.floor(time.time() - self.startTime)
            if secondsSpent < 3600:
                timeStr = time.strftime('%M:%S', time.gmtime(secondsSpent))
            elif secondsSpent < 3600*24:
                timeStr = time.strftime('%H:%M:%S', time.gmtime(secondsSpent))
            else:
                days = secondsSpent // (3600*24)
                seconds = secondsSpent % (3600*24)
                timeStr = int(days) + time.strftime('days %H:%M:%S', time.gmtime(seconds))
            string += ' (' + timeStr + ')'
        string = '\r ' + string + '\r'
        print(string, end='')

    def __del__(self):
        print()

def growBar(bar:ProgressBar, *args, **kw):
    def decorator(func):
        def wrapper(*args1, **kw1):
            res = func(*args1, **kw1)
            bar.grow(*args, **kw)
        return wrapper
    return decorator
