function check_and_create_dirs(filepath)
    % Extrai o caminho do diretório do nome completo do arquivo
    [directory, ~, ~] = fileparts(filepath);
    
    % Verifica se o diretório existe
    if ~exist(directory, 'dir')
        % Se o diretório não existir, cria-o
        mkdir(directory);
        fprintf('Directory created: %s\n', directory);
    end
end