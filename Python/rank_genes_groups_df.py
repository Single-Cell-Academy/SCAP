def rank_genes_groups_df(adata, key='rank_genes_groups'):
    # create a data frame with columns from .uns['rank_genes_groups'] (eg. names, 
    # logfoldchanges, pvals). 
    # Ideally, the list of columns should be consistent between methods
    # but 'logreg' does not return logfoldchanges for example
    import pandas as pd
    dd = []
    groupby = adata.uns['rank_genes_groups']['params']['groupby']
    for group in adata.obs[groupby].cat.categories:
        cols = []
        # inner loop to make data frame by concatenating the columns per group
        for col in adata.uns[key].keys():
            if col != 'params':
                   cols.append(pd.DataFrame(adata.uns[key][col][group], columns=[col]))
        
        df = pd.concat(cols,axis=1)
        df['group'] = group
        dd.append(df)

    # concatenate the individual group data frames into one long data frame
    rgg = pd.concat(dd)
    rgg['group'] = rgg['group'].astype('category')
    return rgg.set_index('group')