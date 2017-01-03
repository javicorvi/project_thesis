nsus_header = ["RUNS\NSUS","1.0","2.0","3.0","4.0","5.0"]
from prettytable import PrettyTable

def get_string(self, **kwargs):

    """Return string representation of table in current state.
       ...
    """

def generate_table(file_path, file_result_path, beta, run, nsus):
    f=open(file_result_path, 'w')
    for b in beta:
        t = PrettyTable(nsus_header)
        print 'Results For Beta: ' +  b
        for r in run: 
            data=[None] * (len(nsus)+1)
            i=0
            data[i]=r
            for n in nsus:
                source = open(file_path, 'r')
                i=i+1
                for line in source:
                    if (('beta'+b in line) & ('nsus'+n in line) & ('runs'+r+'.' in line)):
                        str,value=line.split(",")
                        data[i]=value
                source.close()        
            t.add_row(data)  
        print t            
        data_t = t.get_string()
        f.write('Results For Beta: ' +  b + '\n')
        f.write(data_t)
        f.write('\n')
    f.flush()
    f.close()