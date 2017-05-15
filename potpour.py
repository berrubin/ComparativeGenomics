#!/usr/bin/env python2



import multiprocessing

 
class Worker(multiprocessing.Process):
 
    def __init__(self, work_queue, result_queue, func):
 
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.func = func
        
 
    def run(self):
        while not self.kill_received:
            # get a task
#            if self.work_queue.qsize() == 0:
            if self.work_queue.empty():
                break
            else:
                job = self.work_queue.get()
            res = self.func(*job)
            self.result_queue.put(res)
#if other shit is broken you may have to restore the function as below.
#I don't know why this happened but the hka test wouldn't multithread
#until i did this change. work_queue.empty() kept returning True when
#it was definitely not empty. Frustrating.
"""
            if self.work_queue.empty():
                break
            else:
                #job = self.work_queue.get_nowait()
                job = self.work_queue.get()
 
            # the actual processing
            print "process"
            res = self.func(*job)
 
            # store the result
            self.result_queue.put(res)

"""
