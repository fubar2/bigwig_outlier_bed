## Find and mark BigWig peaks to a bed file for display

In the spirit of DeepTools, but finding contiguous regions where the bigwig value is either above or below a given centile.
0.99 and 0.01 for example. [pybedtools](https://github.com/jackh726/bigtools) is used for the bigwig interface.

These quantile cut point values are found and applied over each chromosome using some [cunning numpy code](http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html)

![image](https://github.com/fubar2/bigwig_peak_bed/assets/6016266/cdee3a2b-ae31-4282-b744-992c15fb49db)

![image](https://github.com/fubar2/bigwig_peak_bed/assets/6016266/59d1564b-0c34-42a3-b437-44332cf1b2f0)


### Note on quantiles per chromosome rather than quantiles for the whole bigwig

It is just not feasible to hold all contigs in the entire decoded bigwig in RAM to estimate quantiles. It would be
better to sample across all chromosomes so as not to lose any systematic differences between them - the current method will hide those
differences unfortunately. Sampling might be possible.
