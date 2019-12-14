import os, logging

TAG = 'txtparsing - '

# This class works with data in txt files
# Uses a template, can work with variety of files in the same format to extract/filter data
class DataWorker:
    # Pulls information from encoding template
    def __init__(self, encoding_template_path):
        try:
            with open(encoding_template_path) as reader:
                line = reader.readline().rstrip('\n')
                self.template = line.split('|')
                
                line = reader.readline().rstrip('\n')
                self.labels = line.split('|')

                reader.close()
        except:
            print('Exception in __init__')

    # Replace an old file with a new one
    @staticmethod
    def replace(old_file, new_file):
        os.remove(old_file)
        os.rename(new_file, old_file)
        return None

    # Prints lines in the range [start, end] for given file
    @staticmethod
    def print_lines(start, end, input_file):
        with open(input_file, 'r') as reader:
            for i, line in enumerate(reader):
                if i >= start and i <= end:
                    print(line)

    # Saves a range of lines to a different file
    @staticmethod
    def read_save(start, end, input_file, output_file):
        writer = open(output_file, 'w')
        with open(input_file, 'r') as reader:
            for i, line in enumerate(reader):
                if i >= start and i <= end:
                    writer.write(line)
        writer.close()

    # Converts a string w/ '+/-' in front to a float
    @staticmethod
    def str_to_flt(str):
        try:
            if str[0] is '+':
                return float(str[1:])
            else:
                return -1 * float(str[1:])
        except:
            print('Error in converting from string to floating point')
            return None

    # Converts float to a string w/ '+/-' in front
    @staticmethod
    def flt_to_str(flt, length):
        if flt > 0:
            term = '+' + str(flt)
        else:
            term = str(flt)
        
        while len(term) != length:
            term += ' '
        return term

    # Saves a list of lines to a text file
    @staticmethod
    def save_lines(output_file, line_list, range):
        with open(output_file, 'w') as writer:
            for line in line_list:
                writer.write(line[range[0]:range[1]] + '\n')

    # Returns a list of terms in a line
    def parse_line(self, line):
        vals = []
        cursor = 0
        for term in self.template:
            if term is '':
                cursor += 1
            else:
                vals.append(line[cursor:cursor + int(term)])
                cursor = cursor + int(term)
        return vals

    # Read a file and save, applying up to 2 range filters
    # Filters are in the format [str:name of variable, low bound, high bound]
    def read_filter(self, input_file, output_file, filter1, filter2=None, prints=False):
        logging.info(TAG+'reading '+input_file+' with the requested filters')
        reader = open(input_file, 'r')
        f1_in = self.labels.index(filter1[0])
        if filter2 is not None:
            f2_in = self.labels.index(filter2[0])
        with open(output_file, 'w') as writer:
            for line in reader.readlines():
                term_list = self.parse_line(line)
                try:
                    term = float(term_list[f1_in][1:])
                    if term_list[f1_in][0] is '-':
                        term *= -1
                except:
                    term = 0
                if term < filter1[2] and term > filter1[1]:
                    if filter2 is not None and term is not 0:
                        try:
                            term = float(term_list[f2_in][1:])
                            if term_list[f2_in][0] is '-':
                                term *= -1
                        except:
                            term = 0
                        if term < filter2[2] and term > filter2[1]:
                            writer.write(line)
                            if prints:
                                print(line)
                    else:
                        writer.write(line)
                        if prints:
                            print(line)
        reader.close()
    
    # Convert a list to a line
    def lst_to_line(self, lst, float_length):
        line = ''
        index = 0
        for term in self.template:
            if term.isdigit():
                if type(lst[index]) is not str:
                    lst[index] = self.flt_to_str(lst[index], float_length)
                line += lst[index]
                index += 1
            else:
                line += ' '
        return line

    # Quicksort a file based on one of the terms, least to greatest
    def quicksort_lg(self, input_file, output_file, label, str_flt_conv=True):
        logging.info(TAG+'starting quicksort of ' + input_file)
        # Convert sorted element from a string to float
        def convert_elements(datalist, index):
            for elem in datalist:
                elem[index] = self.str_to_flt(elem[index])
            return datalist

        # Partitioning, all values before trail are less than pivot value at lst[end]
        def partition(lst, strt, end, int_indx):
            trail = leader = strt
            while leader < end:
                if lst[leader][int_indx] <= lst[end][int_indx]:
                    lst[trail], lst[leader] = lst[leader], lst[trail]
                    trail += 1
                leader += 1
            lst[trail], lst[end] = lst[end], lst[trail]
            return trail
        
        # Actual recursive quicksort
        def quicksort(lst, strt, end, int_index):
            if strt >= end:
                return
            pivot = partition(lst, strt, end, int_index)
            quicksort(lst, strt, pivot-1, int_index)
            quicksort(lst, pivot+1, end, int_index)

        # List of lists, each for substation
        datalist = self.get_vals(input_file, self.labels)
        index = self.labels.index(label)
        if str_flt_conv:
            datalist = convert_elements(datalist, index)

        # Performing the sort
        quicksort(datalist, 0, len(datalist)-1, index)

        with open(output_file, 'w') as writer:
            for term in datalist:
                writer.write(self.lst_to_line(term, 10) + '\n')

    # Returns a list, each index with a list of requested vals from the requested line
    def get_vals(self, input_file, names):
        indecies = []
        vals = []
        for name in names:
            indecies.append(self.labels.index(name))
        with open(input_file, 'r') as reader:
            line = reader.readline()
            while line:
                parsed_list = self.parse_line(line)
                subvals = []
                for index in indecies:
                    subvals.append(parsed_list[index])
                vals.append(subvals)
                line = reader.readline()
        return vals

    def get_vals_lined(self, line, names):
        vals = []
        indecies =[]
        for name in names:
            indecies.append(self.labels.index(name))

        parsed_list = self.parse_line(line)
        for index in indecies:
            vals.append(parsed_list[index])
        return vals