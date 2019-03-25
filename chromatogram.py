from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt



class Chromaotogram:
    """
    :param abi: abi file to parse
    :param show_range: seqence range to show
    :param rev_complement: bool
    :param figsize: tuple, matplotlib figsize
    :param filename: save figure with filename. Default: None, discard.
    
    All of the data necessary for the traces that are conventionally displayed 
    are in the DATA9 through DATA12 channels. 

        >>>from Bio import SeqIO
        >>>record = SeqIO.read("sanger.seq.ab1",format='abi')
        >>>record.annotations['abif_raw'].keys()
        >>>channel_base_order = record.annotations['abif_raw']['FWO_1']
        
    DATA.9-DATA.12: Vectors containing the signal intensities for each channel.
    FWO_1: A string containing the base corresponding to each channel. 
           For example, if it is "ACGT", then DATA9 = A, DATA10 = C, DATA11 = G and DATA12 = T.
    PLOC2: Peak locations as an index of the trace vectors.
    
    _SPCTAGS = [
                "PBAS2",  # base-called sequence
                "PCON2",  # quality values of base-called sequence
                "SMPL1",  # sample id inputted before sequencing run
                "RUND1",  # run start date
                "RUND2",  # run finish date
                "RUNT1",  # run start time
                "RUNT2",  # run finish time
                # NOTE: The following are used for trace data
                "PLOC2",  # position of peaks
                "DATA1",  # channel1 raw data
                "DATA2",  # channel2 raw data
                "DATA3",  # channel3 raw data
                "DATA4",  # channel4 raw data
                "DATA9",  # channel1 analyzed data
                "DATA10",  # channel2 analyzed data
                "DATA11",  # channel3 analyzed data
                "DATA12",  # channel4 analyzed data
                "FWO_1",  # base order for channels
               ]
    """
    def __init__(self, abi, show_range=(0,50), rev_complement=False, figsize=(10,5), **kwargs):
        record = SeqIO.read(abi, format='abi')
        ## {G: DATA9, A: DATA10, T:DATA11, C: DATA12}
        self.channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
        self.channel_base = list(record.annotations['abif_raw']['FWO_1'])
        self.channel_colors=['black','green','red','blue']
        self.color2base = { k:v for k, v in zip(channel_colors, channel_base)}
        self.base2color = { k:v for v, k in zip(channel_colors, channel_base)}
        self.base_peak_index = record.annotations['abif_raw']['PLOC2']      
        self.trace = defaultdict(list)    
        self.figsize = figsize
        self.start, self.end = show_range

        if rev_complement:
            ## {G: DATA9, A: DATA10, T:DATA11, C: DATA12}
            self.sequence = record.reverse_complement().seq
            self.ch_start, self.ch_end = self.base_peak_index[start], self.base_peak_index[end]
            # CTAG
            self.channel_colors=['blue','red','green','black']
            for c in channels: 
                self.trace[c] = list(reversed(record.annotations['abif_raw'][c]))
            # reversed the peak location index of reversed complement sequence
            ch_len = len(trace[c])
            base_peak_index = list(reversed([ch_len-j-1 for j in self.base_peak_index]))
        else:
            self.sequence = record.seq
            self.ch_start, self.ch_end = self.base_peak_index[self.start], self.base_peak_index[self.end] 
            for c in channels: 
                self.trace[c] = record.annotations['abif_raw'][c]


    def plot(self, filename=None):
        """plotting"""
        fig, ax = plt.subplots(figsize=self.figsize)
        for ch, c in zip(self.channels, self.channel_colors): 
            ax.plot(range(self.ch_start, self.ch_end),self.trace[ch][self.ch_start:self.ch_end], 
                    color=c, lw=1, label=self.color2base[c])
            ax.fill_between(range(ch_start,ch_end), 0, trace[ch][ch_start:ch_end], 
                            facecolor=c, alpha=0.125)
        # Plot bases at peak positions
        for i, peak in zip(range(self.start, self.end), self.base_peak_index[self.start:self.end]):
            ax.text(peak, -25, self.sequence[i],
                    color=self.base2color[sequence[i]],
                    horizontalalignment="center")    
        #ax.set_yticklabels([])
        ax.legend(bbox_to_anchor=(1.02, 0.5), loc=2, borderaxespad=0.)
        
        if filename is not None:
            fig.savefig(filename)
            return


