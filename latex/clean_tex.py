import glob
for file in glob.glob('*.tex'):
    with open(file, 'r', encoding='utf-8') as f:
        content = f.read()

    content = content.replace('\u2010', '-')
    content = content.replace('\u2011', '-')
    content = content.replace('\u2012', '-')
    content = content.replace('\u2013', '--')
    content = content.replace('\u2014', '---')
    content = content.replace('\u2015', '---')
    content = content.replace('\u2018', "'")
    content = content.replace('\u2019', "'")
    content = content.replace('\u201c', '"')
    content = content.replace('\u201d', '"')

    with open(file, 'w', encoding='utf-8') as f:
        f.write(content)
