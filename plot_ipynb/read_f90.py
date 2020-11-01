#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
    
def read_parameters(filename, argname, argtype=float):
    """
    Read the value of argname from filename, in whose form ' argname = value'.
    Note that there is two spaces, one just before 'argname' and the other between 'argname' and '='.

    Parameters
    ----------
    filename : str
    argname : str
    argtype : Numeric types {int, float}

    Returns
    -------
    arg : argtype
        The value of 'argname' in numeric type 'argtype'.
    """
    import re
    with open(filename) as f:                             # ASCIIファイルを開く。
        for line in f.readlines():                        # ファイルの各行を読み込み、
            if not (line.strip()).startswith("!"):        # !で始まる行はFortranコメント行なので除外。
                if line.find(" " + argname + " =") != -1: # " argname ="という文字列を含むかどうか判定。
                    arg = line
    arg=re.sub(r'.+' + argname + ' =','',arg) # "argname ="以前の文字を削除。(正規表現： . 改行以外の任意の文字列, + 直前の文字のくり返し)
    arg=re.sub(r'[,!].+\n','',arg)            # コロン","または感嘆符"!"以降の文字と改行コード\nを削除。（正規表現： [abc] "a"または"b"または"c"）
    arg=re.sub(r'd','e',arg)                  # Fortran の倍精度実数を表す d の文字を Pythonの実数でも使える e に置き換える。
    arg=re.sub(r'_DP','e0',arg)               # Fortran の倍精度実数を表す _DP の文字を Pythonの実数でも使える e0 に置き換える。
    if (argtype==str):
        arg=re.sub(r'["\']','',arg).strip()   # Fortran文字列の""や''を削除。
    else:
        arg=argtype(arg)                      # 文字列型を argtype型に変換。
    return arg



def read_time(filename):
    """
    Read the value of time from filename, in specific form 'time= value'.
    Note that there is no space between 'time' and '='.

    Parameters
    ----------
    filename : str

    Returns
    -------
    time : float
        The value of 'time= value'.
    """
    import re
    with open(filename) as f:            # ASCIIファイルを開く
        for line in f.readlines():       # ファイルの各行を読み込み、
            if line.find("time=") != -1: # "time="という文字列を含むかどうか判定。
                arg = line
    arg=re.sub(r'.+time=','',arg)  # "time="以前の文字を削除。(正規表現： . 改行以外の任意の文字列, + 直前の文字のくり返し)
    arg=re.sub(r'[,!].+\n','',arg) # コロン","または感嘆符"!"以降の文字と改行コード\nを削除。（正規表現： [abc] "a"または"b"または"c"）
    time=float(arg)
    return time



#関数の動作確認。
if __name__ == '__main__':
    nx = read_parameters("../src/gkvp_header.f90", "nx", int)
    global_ny = read_parameters("../src/gkvp_header.f90", "global_ny", int)
    global_nz = read_parameters("../src/gkvp_header.f90", "global_nz", int)
    global_nv = read_parameters("../src/gkvp_header.f90", "global_nv", int)
    global_nm = read_parameters("../src/gkvp_header.f90", "global_nm", int)
    nprocw = read_parameters("../src/gkvp_header.f90", "nprocw", int)
    nprocz = read_parameters("../src/gkvp_header.f90", "nprocz", int)
    nprocv = read_parameters("../src/gkvp_header.f90", "nprocv", int)
    nprocm = read_parameters("../src/gkvp_header.f90", "nprocm", int)
    nprocs = read_parameters("../src/gkvp_header.f90", "nprocs", int)
    print(nx, global_ny, global_nz, global_nv, global_nm, nprocw, nprocz, nprocv, nprocm, nprocs)

#他の Python プログラムからインポートできるようにするためには、ノートブック(.ipynb)から Pythonスクリプト(.py)に変換する必要がある。
#ターミナル上で $ jupyter nbconvert --to python read_f90.ipynb


# In[ ]:




