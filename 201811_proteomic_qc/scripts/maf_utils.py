from collections import namedtuple
from pathlib import Path
import gzip


class MAF:
    """
    General purpose of MAF reader.

    It assumes the file begins wtih some comments lines (starting with `#`),
    then the column header, and the actual variant records in TSV format.

    Arguments:
        pth (pathlib.Path or str): Path object to the MAF file.
    """
    def __init__(self, pth):
        pth = Path(pth)
        self.pth = pth
        if pth.suffix == '.gz':
            self._file = gzip.open(str(pth), 'rt')
        else:
            self._file = open(str(pth))

        # A reader wrapping the underlying file object
        # which also returns the line number
        self._reader = enumerate(self._file, 1)

        # Header comments appear before column
        self.header_comments = []
        self.raw_columns = self.read_header()

        # Set up columns
        self.columns = self.make_columns(self.raw_columns)
        self._record_cls = self.make_record_class()

    def read_header(self):
        """
        Read the header comments and return the parsed the column header.
        """
        line_no, line = next(self._reader)
        while line.startswith('#'):
            self.header_comments.append(line.rstrip('\n'))
            line_no, line = next(self._reader)

        # Treat the first noncomment line as columns
        return line.rstrip('\n').split('\t')

    def make_columns(self, raw_columns):
        """Define the columns a variant record should store."""
        return [c.lower() for c in raw_columns]

    def make_record_class(self):
        """Define the record class."""
        return namedtuple(f'MAFRecord', self.columns)

    def make_record(self, vals):
        """Given the MAF record values from a row, construct the record"""
        return self._record_cls._make(vals)

    def __iter__(self):
        return self

    def __next__(self):
        line_no, line = next(self._reader)
        cols = line.rstrip('\n').split('\t')
        return self.make_record(cols)


class GDCMAF(MAF):
    """
    GDC MAF reader

    Based on the file name, it will determine the cancer type and variant caller
    of the MAF and passes them as additional columns.
    """
    def __init__(self, pth, caller):
        super().__init__(pth)
        self.caller = caller

    def make_columns(self, raw_columns):
        columns = super().make_columns(raw_columns)
        # Add caller information
        columns.append('caller')
        return columns

    def make_record(self, vals):
        # Add caller information
        r = self._record_cls(*vals, self.caller)
        # Replace chromosome
        if r.chromosome == 'chrM':
            chrom = 'MT'
        else:
            chrom = r.chromosome[3:]
        r = r._replace(chromosome=chrom)
        return r