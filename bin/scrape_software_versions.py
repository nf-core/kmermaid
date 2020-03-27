#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

<<<<<<< HEAD
regexes = {
    'nf-core/nf-kmer-similarity': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Sourmash': ['v_sourmash.txt', r"sourmash version (\S+)"],
}
results = OrderedDict()
results['nf-core/nf-kmer-similarity'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Sourmash'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-core/nf-kmer-similarity-software-versions'
section_name: 'nf-core/nf-kmer-similarity Software Versions'
section_href: 'https://github.com/czbiohub/nf-kmer-similarity'
=======
# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/kmermaid': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['nf-core/kmermaid'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/kmermaid Software Versions'
section_href: 'https://github.com/nf-core/kmermaid'
>>>>>>> upstream/TEMPLATE
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
<<<<<<< HEAD
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.txt', 'w') as f:
=======
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
>>>>>>> upstream/TEMPLATE
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
