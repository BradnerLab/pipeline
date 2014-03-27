#!/usr/bin/env python

import platform
import smtplib
import socket
import subprocess

def system_info():

    # we might want to switch to something more comprehensive and robust like 
    # https://support.hyperic.com/display/SIGAR/Home
    # but for now I intend to just manually interpret the output for each unique fqdn
    # (there should only be a few different systems during beta testing perf. tracking)

    info = "platform:\nsystem=%s,release=%s\n\ncpu:\n" % (platform.system(), platform.release())
    try:
        with open ("/proc/cpuinfo", "r") as cpuinfo:
            info += cpuinfo.read()
    except:
        try:
            info += subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string",
                "hw.cpufrequency", "hw.physicalcpu", "hw.logicalcpu"])
        except:
            pass
    info += "\nmemory:\n"
    try:
        info += subprocess.check_output(["free", "-g"])
    except:
        try:
            info += subprocess.check_output(["sysctl", "-n", "hw.memsize"])
        except:
            pass
    
    # todo: try logging file system speed, but this might be difficult and 
    #       complicated by file system ram caching

    return info 

def share(app_name, version, log):
    gmail_user = "bamliquidator.perf@gmail.com"
    bad_variable_name = "supersecret".replace("s", "5")

    # please don't misuse bad_variable_name.  if it is misused, I'll just reset it,
    # stop using it, and setup a secure HTTP server to replace it

    # gmail limits 500 emails per day -- that should be plenty

    # smtp code based on http://stackoverflow.com/a/12424439/1007353

    fqdn = socket.getfqdn()

    recipient = ['dfbradnerlab@gmail.com']
    subject = fqdn + " - %s %s - performance measurement" % (app_name, version)

    body = "summary:\napp_name=%s,version=%s,fqdn=%s\n\nlog:\n%s\n%s" % (
        app_name, version, fqdn, log, system_info())

    # Prepare actual message
    message = """\From: %s\nTo: %s\nSubject: %s\n\n%s
    """ % (gmail_user, ", ".join(recipient), subject, body)

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.ehlo()
        server.starttls()
        server.login(gmail_user, bad_variable_name)
        server.sendmail(gmail_user, recipient, message)
        server.close()
    except Exception as e:
        print "failed to send mail:", e

if __name__ == "__main__":
    share("tracker", 1.0, "sample log")

'''
   The MIT License (MIT) 

   Copyright (c) 2014 John DiMatteo (jdimatteo@gmail.com)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE. 
'''
