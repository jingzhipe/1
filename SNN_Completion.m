function [Res, Par] = SNN_Completion(E, Omega, Par)
    % Simple implementation of SNN (Smooth Nuclear Norm) completion
    lambda = Par.lambda;
    maxIter = Par.Iter;
    Res = E;

    for iter = 1:maxIter
        for k = 1:size(E, 3)
            Ek = Res(:, :, k);
            [U, S, V] = svd(Ek, 'econ');
            S = diag(S);
            S = max(S - lambda, 0);
            Res(:, :, k) = U * diag(S) * V';
        end
    end
end
