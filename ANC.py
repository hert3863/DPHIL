# File         : ANC.py
# Author       : Neil Massey
# Created      : 09/06/04
# Modified     : 09/06/04
# Version      : 0.2
# Purpose      : Alpha numeric counter e.g. count 0000 to 0010 would go through 000a, 000z etc
# Program args : None

###############################################################################

def IncASCII(c):
    # increment an ascii alphanumeric digit so that it lies in 0..9, a..z, A..Z
    i = ord(c)
    i = i + 1
    if i == 58:     # '9' + 1
        i = 97      #   'a'
    elif i == 123:  # 'z' + 1
        i = 48      #   '0'
    elif i == 91:   # 'Z' + 1
        i = 48      #   '0'
    return chr(i)

###############################################################################

def IncString(s, pos):
    # increment all the ascii characters in the string
    c = IncASCII(s[pos])
    front = s[:pos]
    back  = s[pos+1:]
    all = front + c + back
    if (c == "0") and (not pos == 0):
        all = IncString(all, pos-1)
    return all

###############################################################################

def EnumerateASCII(a):
    # convert a character to sensible digit i.e. make 0=0, 9=9, a=10, z=35
    d = ord(a)
    if d > 96:          # 'a'
        d = d - 87
    elif d > 64:        # 'A'
        d = d - 55
    elif d > 47:        # '0'
        d = d - 48

    return d

###############################################################################

def NumberBetween(s1, s2):
    # calculate the number of ascii characters between s1 and s2
    # i.e. the number of times next has to be called for s1 == s2

    # check the strings are the same length
    l = len(s1)
    if not l == len(s2):
        return 0

    # loop thru and calculate s2[i] - s1[i] * (l - 1 - i) * 36
    d = 0
    for i in range(0, l):
        o2 = EnumerateASCII(s2[i])
        o1 = EnumerateASCII(s1[i])
        t = (o2 - o1)
        t = t * 36 ** (l - i - 1)
        d = d + t
    return d

def IntValue(s):
    d = 0
    l = len(s)
    for i in range(0, l):
        o = EnumerateASCII(s[i])
        t = o * 36 ** (l - i - 1)
        d = d + t
    return d

###############################################################################

class ANC:
    def __init__(self, width = 4):
        self.m_width = width

    def Start(self, start_number = "0"):
        if len(start_number) < self.m_width:
            start_number = start_number.zfill(self.m_width)
        self.m_current_value = start_number

    def Next(self):
        # get the last digit in ascii form
        self.m_current_value = IncString(self.m_current_value, self.m_width-1)

    def Get(self):
        return self.m_current_value

    m_current_value = "0"
    m_width = 4
