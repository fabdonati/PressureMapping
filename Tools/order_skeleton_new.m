function br = order_skeleton_new( sk, res )

% sk = skeleton
P = sk(1,:);
sk(1,:) = [];

[br1, sk] = get_branch( [P;sk], res );
 br2      = get_branch( [P;sk], res );

br1 = [br1(end:-1:1,:); P];

br = [br1; br2];

    function [br, sk] = get_branch(sk,res)
        br = [];
        dist(1) = 1000000;
        while size(sk,1) > 1
            for j = 2 : size(sk,1)
                dist(j) = norm(sk(1,:)-sk(j,:));
            end
            ind = find(dist<res);
            while length(ind) > 1
                [~,i] = max(dist(ind));
                ind(i) = [];
            end
            if isempty(ind)
                sk(1,:) = [];
                break
            end
            dist(2:end) = [];
            br = [br; mean(sk(ind,:),1)];
            sk(1,:) = mean(sk(ind,:),1);
            sk(ind,:) = [];
        end
        sk(1,:) = [];
    end
end