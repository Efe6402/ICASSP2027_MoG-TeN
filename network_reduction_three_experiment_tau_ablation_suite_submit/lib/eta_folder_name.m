function tag = eta_folder_name(eta)

    tag = ...
        sprintf('eta_%0.6f',eta);

    tag = ...
        strrep(tag,'.','p');
end
