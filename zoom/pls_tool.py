# -*- coding: utf-8 -*-
"""
Functionality for performing cross-validation PLS-R and
traditional imaging-transcriptomics paradigm.
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import norm
import sklearn
from tqdm import tqdm
from sklearn import model_selection, linear_model
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error
from joblib import Parallel, delayed
from sklearn import model_selection, linear_model
import multiprocessing
from typing import Union, Tuple

def optimal_component_eval(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    ncomps: Union[np.ndarray, list],
    cv: int, seed: int
) -> Tuple[list, list, pd.DataFrame]:
    
    """
    Evaluate optimal component for PLS-R through cross-validation strategy.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    ncomps: np.ndarray or list
        Optimal component number candidates.
    cv : int
        Number of folds in cross-validation.
    seed : int
        Random seed to control the dataset split.

    Returns
    -------
    best_comps & r : list
        Local optimal component numbers during the current dataset splition
        and the corresponding Pearson's correlation between the predicted SBP
        and actual SBP.
    preds : pd.DataFrame
        Predictions of cross-validation PLS-R.
        
    References
    ----------
    Wang, Y. et al. Spatio-molecular profiles shape the human 
    cerebellar hierarchy along the sensorimotor-association axis. 
    Cell Rep. 43, 113770 (2024).
    """
    
    # Initialization
    cv_inner = cv
    best_comps = []
    r = []
    preds = pd.DataFrame(np.zeros((len(SBP))))
    preds.index = expression.index
    
    # Split samples with KFold split strategy
    sel = model_selection.KFold(n_splits=cv, shuffle=True, random_state=seed)
    sel_inner = model_selection.KFold(n_splits=cv_inner, shuffle=True, random_state=seed)
    
    for tr_ind, te_ind in sel.split(expression):
        mse = []
        expression_cur = expression.iloc[tr_ind,:]
        for tr_i, te_i in sel_inner.split(expression_cur):
            mse.append([])
            # Inner traing sets
            tr_rows = expression.index[tr_ind][tr_i]
            tr_set = pd.DataFrame(expression.loc[tr_rows,:], copy=True)
            tr_y = pd.DataFrame(SBP).loc[tr_rows,:]
            # Inner test sets
            te_rows = expression.index[tr_ind][te_i]
            te_set = pd.DataFrame(expression.loc[te_rows,:], copy=True)
            te_y = pd.DataFrame(SBP).loc[te_rows,:]
            
            # For all optimal component candidate, evaluate their performance
            for nc in ncomps:
                clf = PLSRegression(n_components=nc)
                pred_y = clf.fit(tr_set,tr_y).predict(te_set).flatten()
                mse[-1].append(mean_squared_error(te_y,pred_y))
                
        # Pick current best component number 
        mse = np.sum(mse,axis=0)
        nc = ncomps[np.argmin(mse)]
        best_comps.append(nc)
        #refit full model
        clf = PLSRegression(n_components=nc)
        # Training sets
        tr_rows = expression.index[tr_ind]
        tr_set = pd.DataFrame(expression.loc[tr_rows,:], copy=True)
        tr_y = pd.DataFrame(SBP).iloc[tr_ind,:]
        # Test sets
        te_rows = expression.index[te_ind]
        te_set = pd.DataFrame(expression.loc[te_rows,:], copy=True)
        te_y = pd.DataFrame(SBP).iloc[te_ind,:]
        mod = clf.fit(tr_set,tr_y)
        preds.loc[te_rows,:] = mod.predict(te_set)
        r.append(stats.pearsonr(mod.predict(te_set).flatten(), te_y.values.flatten())[0])        
    return best_comps,r,preds

def wrapper(args):
    return optimal_component_eval(*args)

def run_component_eval(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    ncomps: Union[np.ndarray, list],
    cv: int, seed: int,
    repeats: int
) -> int:
    
    """
    Evaluate optimal component for PLS-R through cross-validation strategy.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    ncomps: np.ndarray or list
        Optimal component number candidates.
    cv : int
        Number of folds in cross-validation.
    seed : int
        Random seed to control the dataset split.
    repeats : int
        How many times should optimal component evaluation be run?

    Returns
    -------
    best_comp : int
        The optimal component number for current PLS-R model.
        
    References
    ----------
    Wang, Y. et al. Spatio-molecular profiles shape the human 
    cerebellar hierarchy along the sensorimotor-association axis. 
    Cell Rep. 43, 113770 (2024).
    """
    
    # Initialization
    cores = multiprocessing.cpu_count()
    #pool = multiprocessing.Pool(processes=cores)
    tasks = []
    results = []
    
    # Evaluate optimal component number with parallel computation
    for x in range(1, repeats+1):
        current_seed = seed*x
        task = (expression, SBP, ncomps, cv, current_seed)
        tasks.append(task)        
    #ans = pool.starmap(optimal_component_eval,tasks)
    with multiprocessing.Pool(processes=cores) as pool:
        ans = list(tqdm(pool.imap_unordered(wrapper, tasks), total=len(tasks)))
    for n in range(len(ans)):
        results.append(pd.DataFrame({
                       'Best_comp':ans[n][0],
                       'Correlation':ans[n][1]}))
    
    # Pick optimal component number
    comp_stat = pd.concat(results, ignore_index=True)
    best_comp = comp_stat['Best_comp'].value_counts().index[0]
    return best_comp

def model_eval(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    best_comp: int, 
    cv: int, repeats: int, 
    seed: int, n_jobs: int
) -> Tuple[pd.DataFrame, list]:
    
    """
    Evaluate the final model performance for the PLS-R model.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    best_comp : int
        The optimal component number for current PLS-R model.
    cv : int
        Number of folds.
    repeats : int
        How many times should model performance evaluation be run?
    seed : int
        Random seed to control the dataset split.
    n_jobs : int
        Number of cores used for parallel computation.

    Returns
    -------
    preds_rep : pd.DataFrame
        The predicted SBP values across repeats.
    scores : list
        Model performances across repeats, as measured by Pearson's correlation.
        
    References
    ----------
    Wang, Y. et al. Spatio-molecular profiles shape the human 
    cerebellar hierarchy along the sensorimotor-association axis. 
    Cell Rep. 43, 113770 (2024).
    """

    clf = PLSRegression(best_comp)
    
    # Help function to compute model performance for each repeat
    def model_eval_i(i, expression, SBP, clf, cv, seed):
        # Initialization
        preds = np.zeros((len(SBP),))
        sel = model_selection.KFold(n_splits=cv, shuffle=True, random_state=(seed*(i+1)))
    
        for fold, (tr_ind, te_ind) in enumerate(sel.split(expression)):
            # Training sets
            tr_set = expression.iloc[tr_ind]
            tr_y = SBP.iloc[tr_ind]
            # Test sets
            te_set = expression.iloc[te_ind]
            te_y = SBP.iloc[te_ind]
            predictions = clf.fit(tr_set, tr_y).predict(te_set).flatten()
            preds[te_ind] = predictions
        return preds
    
    # Running the parallel processing
    results = Parallel(n_jobs=n_jobs)(
        delayed(model_eval_i)(i, expression, SBP, clf, cv, seed) 
        for i in range(repeats)
    )
    # Extract predictions  from results
    preds_rep = pd.DataFrame([result for result in results]).T
    preds_rep.columns = range(repeats)
    preds_rep.index = expression.index

    scores = []
    for i in range(repeats):
        scores.append(stats.pearsonr(preds_rep.iloc[:,i].values, SBP.values[:,0])[0])
    return preds_rep, scores

def pls_perm(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    SBP_perm: pd.DataFrame, 
    best_comp: int, 
    scores: list, 
    cv: int, 
    seed: int, 
    n_jobs: int
) -> Tuple[float, pd.DataFrame]:
    
    """
    Perform spatial permutation test for PLS-R.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    SBP_perm : pd.DataFrame
        Spatial autocorrelation-preserved null model for SBP.
    best_comp : int
        The optimal component number for current PLS-R model.
    scores : list
        Model performances across repeats, as measured by Pearson's correlation.
    cv : int
        Number of folds.
    seed : int
        Random seed to control the dataset split.
    n_jobs : int
        Number of cores used for parallel computation.

    Returns
    -------
    psa : float
        P-value accounting for spatial autocorrelation.
    scores_perm : list
        Model performances across permutation tests, as measured by Pearson's correlation.
        
    References
    ----------
    Wang, Y. et al. Spatio-molecular profiles shape the human 
    cerebellar hierarchy along the sensorimotor-association axis. 
    Cell Rep. 43, 113770 (2024).
    """

    clf = PLSRegression(best_comp)
    
    # Help function to perform permutation test for each permutation
    def pls_perm_p(p, SBP_perm, expression, SBP, cv, clf, seed): 
        # Initialization
        y = SBP_perm.iloc[:, p]    
        preds_p = pd.DataFrame(np.zeros((len(y),1)), index=expression.index)
    
        sel = model_selection.KFold(n_splits=cv, shuffle=True, random_state=seed)
        for tr_ind, te_ind in sel.split(expression):
            # Training sets
            tr_set, tr_y = expression.iloc[tr_ind], y.iloc[tr_ind]
            # Test sets
            te_set, te_y = expression.iloc[te_ind], y.iloc[te_ind]
            mod = clf.fit(tr_set, tr_y)
            preds_te = mod.predict(te_set)
            preds_p.iloc[te_ind, 0] = preds_te.ravel()
        score_p = stats.pearsonr(preds_p.iloc[:,0], SBP.values[:,0])[0]
        return score_p
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(pls_perm_p)(p, SBP_perm, expression, SBP, cv, clf, seed)
        for p in range(SBP_perm.shape[1])
    )
    scores_perm = results
    
    median_score = np.median(scores)
    p_perm = (scores_perm>median_score).sum()/len(scores_perm)    
    return p_perm, scores_perm

def get_R2(
    mod: sklearn.cross_decomposition._pls.PLSRegression, 
    y: np.ndarray
) -> np.ndarray:
    
    """
    Compute variance explained for each PLS-R component.

    Parameters
    ----------
    mod : sklearn.cross_decomposition._pls.PLSRegression
        PLS-R model fitted to X and y.
    y : np.ndarray
        Dependent variable of PLS-R model.

    Returns
    -------
    R2 : np.ndarray with shape (component,)
        Variance explained for each PLS-R component.
    """
    
    component = mod.n_components
    R2 = np.zeros(component)
    T = mod.x_scores_
    for h in range(component):
        th = np.atleast_2d(T[:,h]).T
        R2h = (th.T@y@y.T@th)/(y.T@y@th.T@th)
        R2[h] = R2h[0][0]
    return R2

def get_vip2(
    mod: sklearn.cross_decomposition._pls.PLSRegression,
    X: np.ndarray, y: np.ndarray
) -> np.ndarray:
    
    """
    Compute variance explained for each PLS-R component.

    Parameters
    ----------
    mod : sklearn.cross_decomposition._pls.PLSRegression
        PLS-R model fitted to X and y.
    X & y : np.ndarray
        Independent and dependent variables of PLS-R model.

    Returns
    -------
    vip : np.ndarray with shape (X.shape[1],)
        Variable importance in projection index for independent variables.
        
    References
    ----------
    Mahieu, B., Qannari, E. M. & Jaillais, B. Extension and 
    significance testing of Variable Importance in Projection (VIP) 
    indices in Partial Least Squares regression and Principal 
    Components Analysis. Chemom. Intell. Lab. Syst. 242, 104986 (2023).
    """
    
    R2 = get_R2(mod, y)
    W = mod.x_weights_
    norm = X.shape[1]/sum(R2)
    vip = norm*np.dot(W**2, R2)
    return vip

def vip_perm(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    SBP_perm: pd.DataFrame, 
    best_comp: int, 
    one_sided: bool,
    n_jobs: int
) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    
    """
    Compute gene-level statistics and corresponding permutation p-values.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    SBP_perm : pd.DataFrame
        Spatial autocorrelation-preserved null model for SBP.
    best_comp : int
        The optimal component number for current PLS-R model.
    one_sided : bool
        If True, infer statistical significance via one-sided
        p-values. Else, use two-sided p-values.
    n_jobs : int
        Number of cores used for parallel computation.

    Returns
    -------
    gene_rep : pd.DataFrame
        Gene-level PLS statistics and permutation p-values.
    RC_perm : pd.DataFrame
        Gene-level regression coefficients (RCs) in permutation tests.
    VIP_perm : pd.DataFrame
        Gene-level Variable importance in projection (VIP) in permutation tests.
        
    References
    ----------
    [1] Wang, Y. et al. Spatio-molecular profiles shape the human 
        cerebellar hierarchy along the sensorimotor-association axis. 
        Cell Rep. 43, 113770 (2024).
    [2] Mahieu, B., Qannari, E. M. & Jaillais, B. Extension and 
        significance testing of Variable Importance in Projection (VIP) 
        indices in Partial Least Squares regression and Principal 
        Components Analysis. Chemom. Intell. Lab. Syst. 242, 104986 (2023).
    """
    
    n_perm = SBP_perm.shape[1]
    # Fit the PLS-R model
    X = expression.values
    y = SBP.values
    clf = PLSRegression(n_components=best_comp)
    mod = clf.fit(X, y)
    
    # Get `RC` and `VIP`
    f_rc = mod.coef_.flatten()
    vip2 = get_vip2(mod, X, y)
    
    # Construct result DataFrame
    gene_rep = pd.DataFrame({'RC': f_rc,
                             'VIP': vip2}, index=expression.columns)
    
    # Help function to perform permutation test for each permutation
    def vip_perm_p(p, X, SBP_perm, best_comp):
        # Initialize PLS-R model
        clf_p = PLSRegression(n_components=best_comp)
        y_p = np.atleast_2d(SBP_perm.iloc[:,p].values).T
        # Fit the permutation null model
        mod_p = clf_p.fit(X, y_p)
        RC_p = mod_p.coef_.flatten()
        VIP2_p = get_vip2(mod_p, X, y_p)
        return RC_p, VIP2_p
    
    # Permutated beta and vip
    results = Parallel(n_jobs=n_jobs)(
        delayed(vip_perm_p)(p, X, SBP_perm, best_comp)
        for p in range(n_perm)
    )
    # Unzip parallel permutation results
    RC_perm = pd.concat([pd.Series(res[0]) for res in results], axis=1)
    RC_perm.index = expression.columns
    VIP_perm = pd.concat([pd.Series(res[1]) for res in results], axis=1)
    VIP_perm.index = expression.columns
    
    # Compute permutation p-values
    if one_sided:
        for g in gene_rep.index:
            RC_g = gene_rep.loc[g,'RC']
            VIP_g = gene_rep.loc[g,'VIP']
            RC_perm_g = RC_perm.loc[g,:]
            VIP_perm_g = VIP_perm.loc[g,:]
            if RC_g > 0:
                condition = RC_perm_g>0
                RC_perm_g = RC_perm_g[condition]
                VIP_perm_g = VIP_perm_g[condition]
                count_rc = sum(RC_perm_g>RC_g)
                count_vip = sum(VIP_perm_g>VIP_g)
            elif RC_g < 0:
                condition = RC_perm_g<0
                RC_perm_g = RC_perm_g[condition]
                VIP_perm_g = VIP_perm_g[condition]
                count_rc = sum(RC_perm_g<RC_g)
                count_vip = sum(VIP_perm_g>VIP_g)                
            gene_rep.loc[g, 'p_perm_rc'] = count_rc/n_perm
            gene_rep.loc[g, 'p_perm_vip'] = count_vip/n_perm        
    else:   
        gene_rep['p_perm_rc'] = [abs(RC_perm.loc[g,:]).ge(abs(gene_rep.loc[g,'RC'])).sum()/n_perm
                                 for g in gene_rep.index]
        gene_rep['p_perm_vip'] = [abs(VIP_perm.loc[g,:]).ge(abs(gene_rep.loc[g,'VIP'])).sum()/n_perm
                                  for g in gene_rep.index]        
    return gene_rep, RC_perm, VIP_perm

def get_Q2(
    SBP: pd.DataFrame,
    preds_rep: pd.DataFrame,
) -> float:
    
    """
    Compute cross-validation R2(Q2).

    Parameters
    ----------
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    preds_rep : pd.DataFrame
        The predicted SBP values across repeats.

    Returns
    -------
    Q2 : float
        The cross-validation R2(Q2).
        
    References
    ----------
    [1] Wang, Y. et al. Spatio-molecular profiles shape the human 
        cerebellar hierarchy along the sensorimotor-association axis. 
        Cell Rep. 43, 113770 (2024).
    [2] Mahieu, B., Qannari, E. M. & Jaillais, B. Extension and 
        significance testing of Variable Importance in Projection (VIP) 
        indices in Partial Least Squares regression and Principal 
        Components Analysis. Chemom. Intell. Lab. Syst. 242, 104986 (2023).
    """
    
    # Compute total sum of squares and residuals
    tss = np.sum((SBP.values-np.mean(SBP))**2)
    residual = preds_rep-SBP.values
    
    # Compute Q2 across repeats
    repeats = preds_rep.shape[1]
    q2_list = []
    for i in range(repeats):
        press = np.sum((residual[i])**2)
        q2 = 1-press/tss
        q2_list.append(q2)
    Q2 = np.median(q2_list)
    return Q2

def boot_pls1(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame,
    best_comp: int,
    n_boot: int,
    seed: int,
    n_jobs: int
):
    """
    Estimate Z-scored PLS1 gene weights with bootstrap strategy.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    best_comp : int
        The optimal component number for current PLS-R model.
    n_boot : int
        Number of bootstrap iteration.
    seed : int
        Random seed to control the resampling.
    n_jobs : int
        Number of cores used for parallel computation.    

    Returns
    pls1_zscore : np.ndarray
        Z-scored PLS1 gene weights estimated through bootstrap.
        
    References
    ----------
    Whitaker, K. J., Vértes, P. E., Romero-Garcia, R & Bullmore, E. T.
    Adolescence is associated with genomically patterned consolidation 
    of the hubs of the human brain connectome. Proc. Natl Acad. Sci. 
    USA 113, 9105-9110 (2016).
    """
    
    n_samples = expression.shape[0]
    # Fit the PLS-R model
    X = expression.values
    y = SBP.values
    clf = PLSRegression(n_components=best_comp)
    mod = clf.fit(X, y)
    # Extract the original weights
    orig_weights = mod.x_weights_[:,0]
    
    # Help function to perform bootstrap PLS1 weights estimation
    def boot_pls1_b(b, X, y, best_comp, n_samples, 
                    orig_weights, seed):
        rng = np.random.RandomState(seed*b)
        # sample indices with replacement
        idx = rng.randint(0, n_samples, size=n_samples)
        X_b = X[idx]
        y_b = y[idx]
        # Fit PLS on bootstrap sample with same n_components
        clf_b = PLSRegression(n_components=best_comp)
        mod_b = clf_b.fit(X_b, y_b)
        boot_weights_b = mod_b.x_weights_[:,0]
        # The sign of PLS components is arbitrary, make sure this aligns between runs
        corr, _ = stats.pearsonr(orig_weights,boot_weights_b)
        if corr < 0:
            boot_weights_b *= -1 
        return boot_weights_b
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(boot_pls1_b)(b,X,y,best_comp,n_samples,orig_weights,seed)
        for b in range(n_boot)
    )    
    boot_weights = np.vstack(results)
    # Bootstrap standard error
    se = boot_weights.std(axis=0, ddof=1)
    pls1_zscore = np.divide(orig_weights, se, out=np.full_like(orig_weights,np.nan), where=se!=0)
    return pls1_zscore

def pls1_perm(
    expression: pd.DataFrame, 
    SBP: pd.DataFrame, 
    SBP_perm: pd.DataFrame, 
    best_comp: int,
    n_boot: int,
    one_sided: bool,
    seed: int,
    n_jobs: int
) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    
    """
    Compute gene-level statistics and corresponding permutation p-values.

    Parameters
    ----------
    expression : pd.DataFrame
        AHBA gene expression matrix.
    SBP : pd.DataFrame
        One-dimensional data frame of spatial brain phenotype.
    SBP_perm : pd.DataFrame
        Spatial autocorrelation-preserved null model for SBP.
    best_comp : int
        The optimal component number for current PLS-R model.
    n_boot : int
        Number of bootstrap iteration.
    one_sided : bool
        If True, infer statistical significance via one-sided
        p-values. Else, use two-sided p-values.
    seed : int
        Random seed to control the resampling.
    n_jobs : int
        Number of cores used for parallel computation.

    Returns
    -------
    gene_rep : pd.DataFrame
        Gene-level Z-scored PLS1 weights and permutation p-values.
    PLS1_perm : pd.DataFrame
        Gene-level Z-scored PLS1 weights in permutation tests.
    """
    
    # Initialization
    n_gene = expression.shape[1]
    n_perm = SBP_perm.shape[1]
    PLS1_perm = np.zeros((n_gene, n_perm))
    
    # Calculate the PLS1 weights for given SBP
    pls1 = boot_pls1(expression, SBP, best_comp, n_boot, seed, n_jobs)
    gene_rep = pd.DataFrame(pls1, index=expression.columns, columns=["PLS1"])
    
    # Calculate the PLS1 weights for each permutation
    for p in tqdm(range(n_perm)):
        y_p = np.atleast_2d(SBP_perm.iloc[:,p].values).T
        y_p = pd.DataFrame(y_p, index=SBP_perm.index)
        pls1_perm_p = boot_pls1(expression, y_p, best_comp,
                                n_boot, seed, n_jobs)
        PLS1_perm[:,p] = pls1_perm_p
    PLS1_perm = pd.DataFrame(PLS1_perm, index=expression.columns)
    
    if one_sided:
        for g in gene_rep.index:
            PLS1_g = gene_rep.loc[g,'PLS1']
            PLS1_perm_g = PLS1_perm.loc[g,:]
            if PLS1_g > 0:
                condition = PLS1_perm_g>0
                PLS1_perm_g = PLS1_perm_g[condition]
                count_pls1 = sum(PLS1_perm_g>PLS1_g)
            elif PLS1_g < 0:
                condition = PLS1_perm_g<0
                PLS1_perm_g = PLS1_perm_g[condition]
                count_pls1 = sum(PLS1_perm_g<PLS1_g)
            gene_rep.loc[g, 'p_perm'] = count_pls1/n_perm
    else:
        gene_rep['p_perm'] = [abs(PLS1_perm.loc[g,:]).ge(abs(gene_rep.loc[g,'PLS1'])).sum()/n_perm
                              for g in gene_rep.index]    
    return gene_rep, PLS1_perm