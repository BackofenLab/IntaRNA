from pathlib import Path
import re
import subprocess
p=Path('/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-246/src/IntaRNA/Matrix.h')
s=subprocess.check_output(["git","show","HEAD:src/IntaRNA/Matrix.h"],cwd=p.parents[2],text=True)
common={
    'size1': 'Number of logical rows.\n@return logical row count',
    'size2': 'Number of logical columns.\n@return logical column count',
    'storageSize': 'Number of allocated elements, including any band padding.\n@return stored element count',
    'clear': 'Reset stored cells to default values without changing the shape.',
    'swap': 'Exchange shape and owned storage.\n@param other matrix to exchange with',
    'operator=': 'Take ownership of another matrix; self-move leaves this matrix unchanged.\n@param other matrix whose storage is moved; emptied on a non-self move\n@return this matrix',
}
for cls in ('Matrix','UpperTriangularMatrix','UpperBandedMatrix'):
    begin=s.index('template<class T>\nclass '+cls+' {')
    end=s.index('\n};',begin)+3
    block=s[begin:end]
    pattern=r'^\t(?:static )?(?:std::size_t |const T &|T &|void |'+cls+r' &)?(?:'+cls+r'|operator\(\)|operator=|size1|size2|storageSize|clear|swap|resize|offset|count|width)\('
    definitions=[]
    edits=[]
    for match in re.finditer(pattern,block,re.M):
        start=match.start()
        # For operator(), the first parentheses belong to the operator name.
        open_paren=block.index('(',start)
        if 'operator()' in match.group(): open_paren=block.index('(',open_paren+2)
        depth=1
        i=open_paren+1
        while depth:
            depth += (block[i]=='(')-(block[i]==')')
            i+=1
        tail=block[i:]
        if tail.lstrip().startswith('= default;'): continue
        body_start=block.index('{',i)
        # A semicolon before the next body indicates a declaration/defaulted method.
        if ';' in block[i:body_start]: continue
        depth=1
        j=body_start+1
        while depth:
            depth += (block[j]=='{')-(block[j]=='}')
            j+=1
        head=block[start:body_start].strip()
        body="\n".join(line[1:] if line.startswith("\t") else line for line in block[body_start:j].split("\n"))
        constructor=head.startswith(cls+'(')
        declaration=re.split(r'(?<!:):(?!:)',head,maxsplit=1)[0].rstrip() if constructor else head
        name=cls if constructor else re.search(r'(operator\(\)|operator=|size1|size2|storageSize|clear|swap|resize|offset|count|width)\(',head).group(1)
        comment=common.get(name)
        if constructor:
            if '&&other' in head:
                comment='Take ownership of another matrix and leave it empty.\n@param other matrix whose storage is moved'
            elif cls=='Matrix':
                comment='Allocate a dense matrix initialized to one value.\n@param rows logical row count\n@param columns logical column count\n@param value initial value for each cell\n@throws std::length_error if the storage size overflows'
            elif cls=='UpperTriangularMatrix':
                comment='Allocate a default-initialized square upper triangle.\n@param rows logical row count\n@param columns logical column count; must equal rows\n@throws std::invalid_argument if the shape is not square\n@throws std::length_error if the storage size overflows'
            else:
                comment='Allocate a default-initialized upper band.\n@param rows logical row count\n@param columns logical column count\n@param lower lower bandwidth; must be zero\n@param upper number of stored superdiagonals\n@throws std::invalid_argument if lower is not zero\n@throws std::length_error if the storage size overflows'
        elif name=='operator()':
            const=' const' in head
            constraint={'Matrix':'within the logical shape','UpperTriangularMatrix':'with i <= j < size2()','UpperBandedMatrix':'within the stored upper band'}[cls]
            comment=('Read a logical cell; structural zeros are returned outside stored cells.' if const and cls!='Matrix' else 'Access a cell '+constraint+'.')
            comment+='\n@param i zero-based row within the logical shape\n@param j zero-based column within the logical shape\n@return '+('read-only' if const else 'mutable')+' reference to the cell'
        elif name=='resize':
            names={'Matrix':('newRows','newColumns'),'UpperTriangularMatrix':('rows','columns'),'UpperBandedMatrix':('rows','newColumns')}[cls]
            comment='Resize, preserving overlapping logical cells when requested.\n@param '+names[0]+' new row count\n@param '+names[1]+' new column count'
            if cls=='UpperTriangularMatrix': comment+='; must equal rows'
            if cls=='UpperBandedMatrix': comment+='\n@param lower lower bandwidth; must be zero\n@param upper number of stored superdiagonals'
            comment+='\n@param preserve retain overlapping cells and default-initialize new cells;\notherwise callers must overwrite cells before reading\n@throws std::length_error if the storage size overflows'
            if cls!='Matrix': comment+='\n@throws std::invalid_argument if the requested shape or band is unsupported'
        elif name=='offset': comment='Offset of a stored upper-triangular cell in the packed allocation.'
        elif name=='count': comment='Validate a square shape and return its packed storage size.'
        elif name=='width': comment='Validate the lower bandwidth and clamp the stored width to the columns.'
        assert comment,(cls,name,head)
        doc='\t/**\n'+''.join('\t * '+line+'\n' for line in comment.splitlines())+'\t */\n'
        edits.append((start,j,doc+'\t'+declaration+';'))
        head="\n".join(line[1:] if line.startswith("\t") else line for line in head.split("\n"))
        head=head.replace('static ','',1) if head.startswith('static ') else head
        head=head.replace('const T &value = T{}','const T &value').replace('bool preserve = true','bool preserve')
        if constructor:
            head=head.replace(cls+'(',cls+'<T>::'+cls+'(',1)
        else:
            if head.startswith(cls+' &'): head=head.replace(cls+' &',cls+'<T> &',1)
            head=head.replace(name+'(',cls+'<T>::'+name+'(',1)
        definitions.append('template<class T>\ninline\n'+head+'\n'+body+'\n')
    assert len(definitions)==({'Matrix':11,'UpperTriangularMatrix':13,'UpperBandedMatrix':12}[cls]),(cls,len(definitions))
    for a,b,replacement in reversed(edits): block=block[:a]+replacement+block[b:]
    block=block.replace('\t'+cls+'() = default;', '\t/** Construct an empty matrix. */\n\t'+cls+'() = default;')
    block=block.replace('\t'+cls+'(const '+cls+' &) = default;', '\t/** Copy the shape and values into independent storage. */\n\t'+cls+'(const '+cls+' &) = default;')
    block=block.replace('\t'+cls+' &operator=(const '+cls+' &) = default;', '\t/** Copy the shape and values into independent storage.\n\t * @return this matrix\n\t */\n\t'+cls+' &operator=(const '+cls+' &) = default;')
    s=s[:begin]+block+'\n\n'+'\n'.join(definitions)+s[end:]
s=s.replace('\tstd::vector<T> values;', '\t//! Owned contiguous storage.\n\tstd::vector<T> values;')
s=s.replace('\tstd::size_t rows = 0, columns = 0;', '\t//! Logical dimensions.\n\tstd::size_t rows = 0, columns = 0;')
s=s.replace('\tstd::size_t n = 0;', '\t//! Logical row and column count.\n\tstd::size_t n = 0;')
s=s.replace('\tMatrix<T> band;', '\t//! Physical rows of the upper band, including row-end padding.\n\tMatrix<T> band;')
s=s.replace('\tstd::size_t columns = 0;', '\t//! Logical column count, independent of stored band width.\n\tstd::size_t columns = 0;')
s=s.replace('\tusing value_type = T;', '\t//! Type of one stored cell.\n\tusing value_type = T;')
p.write_text(s)
