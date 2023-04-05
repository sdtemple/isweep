#     # coalescent process
#     j = n - len(Taus)
#     if len(Taus) <= 0:
#         Taus = np.zeros(0)
#         itr = 0
#     else:
#         itr = Taus[-1]
#     if j > 1:
#         sNe = {k:v for k,v in Ne.items() if k >= itr}
#         if len(sNe) > 0:
#             Etas = varying_Ne_coalescent(j, sNe, ploidy, to_tmrca=to_tmrca)
#             Taus = np.concatenate((Taus,Etas))
#         else:
#             if to_tmrca:
#                 lastN = max(Ne.keys())
#                 scale = Ne[lastN]
#                 Etas = basic_coalescent(j) * lastN + itr
#                 Taus = np.concatenate((Taus,Etas))
#         return Taus
#     else:
#         return Taus

# def wright_fisher(n, Ne, ploidy = 2, to_tmrca=True, cut=1):
#     '''Simulate times in Wright Fisher model with varying population size
#     (After last generation in Ne, assume constant size
#     and apply scaled basic_coalescent)
    
#     Parameters
#     ----------
#     n : int
#         Sample size
#     Ne : dict
#         Effective population sizes
#     ploidy : int
#     to_tmrca : bool
#         Go to TMRCA
#     cut : int
#         Limit to stop WF process
    
#     Returns
#     -------
#     NumPy array
#         Interarrival times in generations
#     '''
#     # initalize
#     n = int(float(n))
#     k = n
#     itr = 0
#     gentimes = np.zeros(n-1)
#     Me = deepcopy(Ne)
#     timer = min(Me.keys())
#     lastG = max(Me.keys())
#     lastN = Me[lastG]
#     cuml = 0
#     while k > 1 and timer <= lastG:
#         # record coalescent events
#         size = Me[timer] * ploidy
#         draw = [randint(0, size - 1) for i in range(k)]        
#         table = {}
#         for i in range(len(draw)):
#             try:
#                 table[draw[i]] += 1
#             except KeyError:
#                 table[draw[i]] = 1
#         events = [val for key, val in table.items() if val >= 2]
#         l = sum(events) - len(events)
#         cuml += l
#         timer += 1
#         gentimes[itr:cuml] = timer
#         k -= l
#         itr += l
#         if k <= cut:
#             return gentimes[gentimes > 0]
#     if to_tmrca:
#         if k > 1: # finish with coalescent
#             finish = basic_coalescent(k) * lastN + lastG
#             gentimes[cuml:] = finish
#     return gentimes[gentimes > 0]

#     minG = min(Ne.keys())
#     Ne = {(k-minG):v for k,v in Ne.items()}

# def nonzero_round(val):
#     val=round(val)
#     if val > 0:
#         return val
#     else:
#         return 1

# def np_nonzero_round(vals):
#     vals=np.round(vals)
#     vals[vals==0]=1
#     return vals

# def np_nonzero_floor(vals):
#     vals=np.floor(vals)
#     vals[vals==0]=1
#     return vals
