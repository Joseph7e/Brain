import csv

class Tagged(object):

    def __init__(self, data, tags):
        assert len(data[0]) == len(tags)

        self.data = data
        self.tags = tags

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.data[self.tags.index(attr)]
        elif isinstance(key, int):
            return self.data[key]
        else:
            raise TypeError()

    def __iter__(self):
        for x in self.data:
            yield dict(zip(self.tags, x))

    def get(self, index):
        return dict(zip(self.tags, self.data[index]))

    def slice(self, start, end):
        return Tagged(
            self.data[slice(start, end)],
            self.tags
        )

    def indices(self, new_tags):
        return map(
            lambda x: self.tags.index(x),
            new_tags
        )

    def condense(self, new_tags):
        assert set(new_tags).issubset(set(self.tags))

        return Tagged([
            map(
                lambda x: col[x],
                self.indices(new_tags)
            ) for col in self.data
        ],
        new_tags)

    def integerize(self, new_tags):
        return Tagged(
            zip(*[
                map(
                    lambda x: int(x),
                    col
                )
                if i in self.indices(new_tags)
                else col 
                for i, col in enumerate(
                    zip(*self.data)
                )
            ]),
        self.tags)


class Genome(object):

    def __init__(self, fasta, misa, edge_length=25):
        self.contigs = self.read_fasta(fasta)
        self.satellites = self.read_misa(misa)

        self.edge_length = edge_length

    def memoize(f):
        memo = {}
        def wrapper(self, *args):
            key = (hash(self.contigs), hash(self.satellites)) + args
            if key in memo:
                return memo[key]
            else:
                val = f(self, *args)
                memo[key] = val
                return val
        return wrapper

    # File -> [[id, seq]]
    def read_fasta(self, fasta):
        return Tagged([
                [
                    contig[0],
                    contig[1].replace('\n', '')
                ]
                for contig in map(
                    lambda x: x.split('\n', 1),
                    open(fasta, 'r')
                        .read()
                        .split('>')[1:]
                )
            ],
            ['id', 'seq']
        )

    # File -> [[id, start, end]]
    def read_misa(self, misa):
        return Tagged(
            list(
                csv.reader(
                    open(misa, 'r'),
                    delimiter='\t'
                )
            )[1:],
            ['id', 'ssr_nr', 'ssr_type', 'ssr', 'size', 'start', 'end']
        ).condense(['id', 'ssr', 'start', 'end']).integerize(['start', 'end'])

    # File -> [[id, leading, trailing]]
    @memoize
    def edges(self):
        return Tagged(
            [
                [
                    contig['id'],
                    satellite['ssr'],
                    contig['seq'][satellite['start'] - self.edge_length: satellite['start']],
                    contig['seq'][satellite['end']: satellite['end'] + self.edge_length]
                ]
                for satellite in self.satellites
                for contig in filter(
                    lambda x: satellite['id'] == x['id']
                        and satellite['start'] >= self.edge_length
                        and satellite['end'] + self.edge_length <= len(x['seq']),
                    self.contigs
                )
            ],
            ['id', 'ssr', 'leading', 'trailing']
        )

    # Genome -> Genome -> [[foreign_id, foreign_ssr, id, ssr, leading, trailing]]
    def match(self, genome):
        assert isinstance(genome, Genome)
        
        return Tagged(
            [
                [edge['id'], edge['ssr']] +
                filter(
                    lambda x: edge['leading'] == x['leading']
                    and edge['trailing'] == x['trailing'],
                    self.edges
                )
                for edge in genome.edges()
            ],
            ['foreign_id', 'foreign_ssr', 'id', 'ssr', 'leading', 'trailing']
        )



A02 = Genome('ANN-A02_Contigs.fasta', 'ANN-A02_Contigs.fasta.misa')
B7 = Genome('ANN-B7_Contigs.fasta', 'ANN-B7_Contigs.fasta.misa')
print A02.match(B7)




