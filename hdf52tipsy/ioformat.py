def rcol(filename="", columns=[], intcolumns=[], linestart=0, verbose=True):
    if(filename=="help" or filename=="" or len(columns)==0):
        print "Usage:"
        print "rcol(filename, [], intcolumns=[], linestart=0)"
        print "Example:"
        print "gid, mstar = rcol(fname, [0,2], [0], linestart=1)"
        return
    if(verbose == True):
        print "Reading File: ", filename
    linecount = 0
    f = open(filename, "r")
    if linestart > 0:
        for i in range(linestart):
            f.readline()
    cols = []
    for i in range(len(columns)):
        cols.append([])
    for line in f: # Each line
        col_i = 0
        spt = line.split()
        for c in columns:
            cols[col_i].append(spt[c])
            col_i = col_i + 1
        linecount = linecount + 1
    f.close()
    if(verbose == True):
        print "Formatting Output:"
    intcols = [0] * len(columns)
    for i in intcolumns:
        intcols[i] = 1
    j = 0
    for col in cols:
        if(intcols[j] == 1):
            for i in range(linecount):
                col[i] = int(col[i])
        if(intcols[j] == 0):
            for i in range(linecount):
                col[i] = float(col[i])
        j = j + 1
    if len(cols) > 1:
        return cols
    else:
        return cols[0]

def printline(datalist, num_of_int=0):
    # CHECKED
    s = ""
    for dat in datalist[:-1]:
        s = s + str(dat) + " "
    s = s + str(datalist[-1]) + "\n"
    return s

def splitline(line, num_of_int=0):
    list1 = line.split()
    data = []
    for s in list1[:num_of_int]:
        data.append(int(s))
    for s in list1[num_of_int:]:
        data.append(float(s))
    return data

def readcol(file, num_of_int=0):
    class columns():
        def __init__(self, ncols, num_of_int):
            self.n = int(ncols); self.num_of_int = num_of_int
            self.cols = []
            for i in range(self.n): self.cols.append([])
    dat = splitline(file.readline(), num_of_int)
    cols = columns(len(dat), num_of_int)
    for i in range(cols.n): cols.cols[i].append(dat[i])
    for line in file:
        dat = splitline(line, num_of_int)
        for i in range(cols.n): cols.cols[i].append(dat[i])
    print cols.n, "columns read, num_of_int =", num_of_int
    return cols
