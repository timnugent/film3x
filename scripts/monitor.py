from pymol import cmd
import threading
import time
import os
import sys
 
class pymol_file_monitor (object) :
  def __init__ (self,
      file_name, 
      time_wait=1) : # time in seconds between mtime check
    self.file_name = file_name
    self.time_wait = time_wait
    self.watch = True # this can be toggled elsewhere to stop updating
    self.mtime = 0
    t = threading.Thread(target=self.check_file)
    t.setDaemon(1)
    t.start()
    print "Watching file %s" % file_name
 
  def check_file (self) :
    while (self.watch) :
      if (os.path.exists(self.file_name)) :
        print "checking..."
        mtime = os.path.getmtime(self.file_name)
        if (mtime > self.mtime) :
          self.mtime = mtime
          print "Re-loading %s" % self.file_name
          time.sleep(2)
          cmd.load(self.file_name, state=1)
          cmd.reset();
          cmd.hide(representation="lines",selection="all")
          cmd.show(representation="cartoon",selection="all")
          view = cmd.get_view()
          camera = list(view)
          camera[0] = -1
          camera[1] = 0
          camera[2] = 0
          camera[3] = 0
          camera[4] = 0
          camera[5] = 1
          camera[6] = 0
          camera[7] = 1
          camera[8] = 0
          cmd.set_view(camera)
      time.sleep(self.time_wait)
 
if (__name__ == "pymol") :
  monitor = pymol_file_monitor("test.pdb")
