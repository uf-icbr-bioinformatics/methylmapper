import subprocess

class Creator():
    """A class that represents a connection to the gdcreate program."""
    gdcreate = "gdcreate"
    proc = None
    pin = None
    pout = None

    def __init__(self, pathname="gdcreate"):
        self.gdcreate = pathname
        self.proc = subprocess.Popen(self.gdcreate, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        self.pin = self.proc.stdin
        self.pout = self.proc.stdout

    def send(self, words):
        """Basic function to send a command to gdcreate."""
        for w in words:
            self.pin.write("{}\n".format(w))
        self.pin.flush()

    def sendCommand(self, *words):
        """Top-level function to send a command to gdcreate and read its reply."""
        return self.sendCommandList(words)

    def sendCommandList(self, words):
        self.send(words)
        reply = self.pout.readline()
        if reply == "bad\n":
            sys.stderr.write("Error in command: {}\n".format(words))
            sys.exit(1)
        return reply.rstrip("\n")

    def terminate(self):
        self.send("ZZ")
        self.proc.terminate()

