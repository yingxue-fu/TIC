import sys
import re
import traceback


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        
        if self.count == 1:
            text = "\r    Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block), int(self.progress * 100), self.status)
            sys.stdout.write(text)
            sys.stdout.flush()
        
        elif self.count == int(0.25*self.total):
            text = "\r    Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block), int(self.progress * 100), self.status)
            sys.stdout.write(text)
            sys.stdout.flush()
            
        elif self.count == int(0.5*self.total):
            text = "\r    Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block), int(self.progress * 100), self.status)
            sys.stdout.write(text)
            sys.stdout.flush()
        
        elif self.count == int(0.75*self.total):
            text = "\r    Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block), int(self.progress * 100), self.status)
            sys.stdout.write(text)
            sys.stdout.flush()
        
        elif self.count == self.total:
            text = "\r    Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block), int(self.progress * 100), self.status)
            sys.stdout.write(text)
            sys.stdout.flush()


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith("feature"):
                parameters["feature_files"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files":
                    parameters[key] = [val]
                else:
                    parameters[key] = val

    return parameters


# Context manager that copies stdout and any exceptions to a log file
class Tee(object):
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout

    def __enter__(self):
        sys.stdout = self

    def __exit__(self, exc_type, exc_value, tb):
        sys.stdout = self.stdout
        if exc_type is not None:
            self.file.write(traceback.format_exc())
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()