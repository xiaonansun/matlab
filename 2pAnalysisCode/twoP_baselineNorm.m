function normMeanMatrix = twoP_baselineNorm(m,r)

normMeanMatrix=bsxfun(@rdivide,bsxfun(@minus,m,min(m,[],2)),max(m,[],2)-min(m,[],2));
end