class OverlapResolver:
    def __init__(self, reads):
        self.reads = reads
        self.overlap_buffer = {}
        self.max_acc = 0

    def get_read_content(self, read_id):
        return self.reads[read_id]

    def get_overlap_by_read_ids(self, left_read_id, right_read_id):
        l = self.get_read_content(left_read_id)
        r = self.get_read_content(right_read_id)
        return self.get_overlap_by_read_contents(l, r)

    def compute_overlap(self, left_read, right_read):
        for i in range(len(left_read)):
            l = left_read[i:]
            size = len(l)
            r = right_read[:size]
            if l == r:
                return l, size
        return "", 0
    
    def get_overlap_by_read_contents(self, left_read, right_read):
        if left_read not in self.overlap_buffer or right_read not in self.overlap_buffer[left_read]:
            overlap = self.compute_overlap(left_read, right_read)[1]
            if left_read not in self.overlap_buffer:
                self.overlap_buffer[left_read] = {}
            self.overlap_buffer[left_read][right_read] = overlap
        return self.overlap_buffer[left_read][right_read]

    def count_reads(self):
        return len(self.reads)
