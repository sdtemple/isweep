# packages
import sys
import os
import gzip
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import binom, norm
from math import sqrt, floor, ceil, exp, log
from random import randint
from copy import deepcopy
import networkx as nx

# my modules
from .outgroups import *
from .favoredalleles import *
from .utilities import *
from .inference import *
from .coalescent import *
from .slow import *
