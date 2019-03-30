#! /usr/bin/env python3

import yaml, subprocess, os

# context to temperarily switch to a different
# working directory.
class working_directory(object):
  def __init__(self,path):
    self.old_dir = os.getcwd()
    self.new_dir = path
  def __enter__(self):
    if self.new_dir is not None:
      os.chdir(self.new_dir)
  def __exit__(self,type,value,traceback):
    if self.new_dir is not None:
      os.chdir(self.old_dir)



with open("manifest.yaml",'r') as f:
  packages = yaml.load(f)



for package in packages:
  directory = packages[package].get("directory",package)
  if not os.path.exists:
    for url in packages[package].get('urls',[]):
      subprocess.check_call(["git","clone", url, directory])
  else:
    with working_directory(directory):
      subprocess.check_call(["git","pull", "origin"])

  checkout = packages[package].get("checkout","master")
  with working_directory(directory):
    subprocess.check_call(["git","checkout", checkout])
