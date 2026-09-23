function blocks = make_component_blocks(N,R)

    base_size = floor(N/R);
    remainder = mod(N,R);

    blocks = repmat( ...
        struct( ...
        'start',0, ...
        'stop',0, ...
        'size',0, ...
        'indices',[]), ...
        R,1);

    cursor = 1;

    for r = 1:R

        this_size = ...
            base_size ...
            +(r<=remainder);

        idx = ...
            cursor:(cursor+this_size-1);

        blocks(r).start = idx(1);
        blocks(r).stop = idx(end);
        blocks(r).size = numel(idx);
        blocks(r).indices = idx;

        cursor = cursor+this_size;
    end
end
