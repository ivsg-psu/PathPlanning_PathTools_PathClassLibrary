% This Wrapper class is used to connect to a database and support basic SQL insert, query, update,
% and delete. 
%
% Author: Liming Gao
% Create Date: 2020-03-17
% =======update=======
% 1.
%======== to do list ============
% 1. figure out the usage of default vlues of properties (done, 04/16/2020)
% 2. for bulk insert, extend the values from numeric matrix to cell
% 3. finish function bulkInsertByFile
% 

classdef Database < handle
    
    properties
        
        database_name
        ip_address = '130.203.223.234' %  'localhost'; %Ip address of server host 
        port = 5432;  % port number 
%         username = 'ivsg_db_user'   % user name for the server
%         password = 'ivsg@DB320' % password
username = 'brennan'; % user name for the server
password = 'ivsg@Reber320'; % password

        db_connection  % database connection, this is assigned in  constructor function 
        
        % Use the PostgreSQL Matlab database driver if you setup the
        % connection using databse toolbox
        driver = 'org.postgresql.Driver'; 
        
        srid = 4326;
        
    end
    
    methods
        
        % Constructor for class. Assign values to properties
        % This uses the internal connect function to create a connection to the database.
        function obj = Database(database_name, ip_address, port, username, password)
            
            if nargin == 1
                obj.database_name = database_name;
                obj = Database.connect(obj);
                
            elseif nargin == 5
                obj.database_name = database_name;
                obj.ip_address = ip_address;
                obj.port = port;
                obj.username = username;
                obj.password = password;
                
                obj = Database.connect(obj);
                
            else
               
                error('Database requires 5 inputs: database name, ip address, port number,  database username and password');
                
            end
            
        end
        
        % Generate query of the database to allow SELECT, INSERT, DELETE, and UPDATE queries.
%         function result = query(obj,q)
%             
%             result = exec(obj.db_connection, char(q));
%             
%         end
        
        % The exec function is not recommended.
        % For SQL statements that return data, use the fetch function or the select function instead. 
        % For other SQL statements, use the execute function instead.
        % Import data into MATLAB workspace from execution of SQL statement
        % results = fetch(conn,sqlquery,opts)
        function result = fetch(obj,q)
            
            result = fetch(obj.db_connection,char(q));
            
        end
        % =============================== function select===============================
        % purpose:           query the data from DB
        % intput:            table = target table, format: string, eg:'gps'
        %                    fields = fields of table, format: cell array, eg: {'id','name'}
        %                    where = query condition, format: cell array, eg: {'id<10','date_added > '2020-02-01''}
        %                    orderby = order by, format: string, eg:'id'
        %                    limit = number of rows to query,format: number
        % output:            r = table format of result
        %               result = struct format of result
        function [r, result, column_names] = select(obj, table, fields, where, orderby, limit)
            
            %%%step1: Build the select statement %%%
            
            if  strcmp(fields, 'all')
                fields = Database.allFields(obj,table);  % find all the fields given a DB table
            end
            
            fields = Database.parseFields(fields);  % this function do nothing right now 
            % Join the fields we want to return
            q = ['SELECT ' strjoin(fields, ', ')];
            % Add the table we are querying
            q = [q ' FROM ' table];
            
            % Add the where clause if it is specified
            if exist('where', 'var')  %Check existence of variable, script, function, folder, or class 
                if ~isempty(where)    
                    q = [q ' WHERE '];
                    for i = 1:length(where)
                        if i ~= length(where)
                            q = [q cat(2, where{i}, ' AND ')]; %#ok<AGROW>
                        else
                            q = [q where{i}]; %#ok<AGROW>
                        end
                    end
%                     q = [q ' WHERE ' strjoin(where, ' AND ')];
                end
            end
            
            % Add the order by clause if it is specified
            if exist('orderby', 'var')                
                if ~isempty(orderby)
                    q = [q ' ORDER BY ' orderby]; 
                end
            end
            
            % Add the limit clause if it is specified
            if exist('limit', 'var')
                if limit > 0
                    q = [q ' LIMIT ' num2str(limit)];
                else
                    error('LIMIT must be a value greater than zero.')      
                end 
            end
            
            % Join everything into one query string
            q = strjoin(string(q));
            
            % And query the database
            r = fetch(obj.db_connection,char(q));
            
            if isempty(r)

                error('Database:emptyQueryResults', char(cat(2,'Input query returns no results. Modify your query!\n\n', {' '}, q)));
                
            end
            
            %%% step2: Extract column names to build structure fields %%%
            result = {}; % create a empty cell array 
            column_names = cell(1,length(fields)); %create a cell array of empty matrices
            for i = 1:length(fields)
                
                % We want to grab the specified column name. There are
                % cases where you can rename the output result 'AS' another
                % name. E.g. 'timestamp AS time', so we want to find
                % 'time'. We will split the string with a
                % delimeter of ' AS '.
                column_name_split = strsplit(fields{i}, ' AS ');
                
                % Field name will be the last element
                field_name = column_name_split(end);
                
                % Remove distinct from before field name. E.g.
                % SELECT DISTINCT trip_id ... -> trip_id
                field_name = strrep(field_name, 'DISTINCT ', '');
                
                field_name = char(field_name);
                
                column_names{i} = field_name;
                    
                % Create an empty struct
                result.(field_name) = r(:,i);
                
            end
            
        end
        % =============================== function insert===============================
        % purpose:           insert one row to database 
        % intput:
        % output:            none
        %               
        function insert(obj, table, fields, values)
           
            fields = Database.parseFields(fields);
           
            %%% Build the select statement %%%
            
            % Start the query by adding the table we want to insert into
            q = ['INSERT INTO' table];
            
            % Add the fields
            q = [q '(' strjoin(fields, ', ') ')'];
            
            % Add the values we want to insert
            q = [q 'VALUES (' values ')'];
            
            % Join everything into one query string
            q = strjoin(q);
            
%             exec(obj.db_connection, char(q));
            curs = execute(obj.db_connection, char(q));
            close(curs);
            
        end
        
        function update(obj, table, update_set, where)
           
            %%% Build the select statement %%%
            
            % Start the query by adding the table we want to insert into
            q = ['UPDATE' table];
            
            % Add the fields
            q = [q 'SET'];
            for i = 1:length(update_set)
                if i ~= length(update_set)
                    q = [q cat(2, update_set{i}, ', ')]; %#ok<AGROW>
                else
                    q = [q update_set{i}]; %#ok<AGROW>
                end
            end
            
            % Add the where clause if it is specified
            if exist('where', 'var')
                
                if ~isempty(where)
                
                    q = [q 'WHERE'];
                    
                    for i = 1:length(where)
                        if i ~= length(where)
                            q = [q cat(2,where{i}, ' AND ')];  %#ok<AGROW>
                        else
                            q = [q where{i}];  %#ok<AGROW>
                        end
                    end
                end
                
            end
            
            % Join everything into one query string
            q = strjoin(q);
            
            %exec(obj.db_connection, char(q));
            execute(obj.db_connection, char(q));
%             curs = execute(obj.db_connection, char(q));
%             close(curs);

        end
        
        
        %------- bulkInsert
        % input: 
        % table, table name,string ,eg: 
        % fields, fields, array, eg:
        % values, values to be inserted, numeric matrix, eg: 
        % http://stackoverflow.com/questions/7019831/bulk-batch-update-upsert-in-postgresql
        function bulkInsert(obj, table, fields, values)
          
            % Start the query by adding the table we want to insert into
            q = ['INSERT INTO' string(table)];
            
            % Add the fields
            q = [q '(' strjoin(fields, ', ') ')'];

            q = [q 'VALUES ('];
            
            for i = 1:length(fields)
               q = [q 'UNNEST(ARRAY['];  %#ok<AGROW>
               q = [q char(strjoin(string(num2str(values(:,i),15)),','))];  %#ok<AGROW>
               q = [q '])'];  %#ok<AGROW>
               
               if i ~= length(fields)
                  q = [q ','];   %#ok<AGROW>
               end
            end
            
            q = [q ')'];
            
            q = strjoin(q);

%             exec(obj.db_connection, char(q));
            curs = exec(obj.db_connection, char(q));
            close(curs);
            
        end
        
        % http://stackoverflow.com/questions/7019831/bulk-batch-update-upsert-in-postgresql
        function bulkUpdate(obj, table, fields, values, where)
            
            % Start the query by adding the table we want to insert into
            q = ['UPDATE' string(table)];
            
            % Add the fields
            q = [q 'SET'];
            for i = 1:length(fields)
                q = [q cat(2,fields{i}, '=data_table.', fields{i})];  %#ok<AGROW>
                
                if i ~= length(fields)
                   q = [q ','];   %#ok<AGROW>
                end
            end

            q = [q 'FROM (SELECT'];
            
            for i = 1:length(fields)
               q = [q 'UNNEST(ARRAY['];  %#ok<AGROW>
               q = [q char(strjoin(string(num2str(values(:,i),15)),','))]; %#ok<AGROW>
               q = [q ']) AS' fields{i}];  %#ok<AGROW>
               
               if i ~= length(fields)
                  q = [q ','];   %#ok<AGROW>
               end
            end
            
            q = [q ') AS data_table'];
            
            q = [q cat(2,'WHERE ',string(table),'.id=data_table.id')];
            
            q = strjoin(q);

%             exec(obj.db_connection, char(q));
            curs = execute(obj.db_connection, char(q));
            close(curs);
            
        end
        
        function bulkInsertByFile(obj, table, fields, values)
            
            % Export Data Using Bulk Insert,https://www.mathworks.com/help/database/ug/export-data-using-bulk-insert.html
            
            % insert into MYsql
            A = {100000.00,'KGreen','06/22/2011','Challengers'};
            A = A(ones(10000,1),:);
            
            A = values;
            
            fid = fopen('tmp.txt','wt');
            for i = 1:size(A,1)
                fprintf(fid,'%10.2f \t %s \t %s \t %s \n', ...
                    A{i,1},A{i,2},A{i,3},A{i,4});
            end
            fclose(fid);
            
            % Run the bulk insert functionality. The SQL statement uses the statement LOCAL INFILE for error handling. For details about this statement, consult the MySQL database documentation.
            
            execute(obj.db_connection,['LOAD DATA LOCAL INFILE ' ...
                ' ''C:\\temp\\tmp.txt'' INTO TABLE ' table  '' ...
                'FIELDS TERMINATED BY ''\t'' LINES TERMINATED ' ...
                'BY ''\n'''])
            
            % Confirm the number of variables in BULKTEST.
            
            results = fetch(conn,'SELECT * FROM BULKTEST');
            results.Properties.VariableNames
            
            % Slight change of plan I am now getting a matrix instead of text file. Sometimes the size of the matrix is large over 100K x 116. Here is my current logic.
            % Steps 1 to 3 are done in in batches of 500 rows per iteration. I am still investigating the optimal way of doing this.
            % 1) Convert the matrix into a string using sprintf.
            % 2) strrep(matrix_to_string, 'NaN', 'NULL'); To replace NaN with Null
            % 3). Write to a text file fprintf statement.
            % 4) Read the text file and perform a Bulk Insert
            
        end
        
        
        function updateGISColumn(obj, table)
           
            update_set = {cat(2,'geography=ST_SetSRID(ST_MakePoint(longitude,latitude),', num2str(obj.srid), ')')};
            update(obj, string(table), update_set);
            
        end
        
        % Find Catalogs and Schemas in the database 
        function tables = ShowTables(obj)
            
            public_schema = sqlfind(obj.db_connection,'','Schema','public');
            tables =  public_schema(strcmp(public_schema.Type,'TABLE'),:);
        end
         % Disconnect from the database.
        function disconnect(obj)
           
            close(obj.db_connection)
            
            fprintf('Disconnected\n')
            
        end
        
    end
    
    methods (Static)
        
        % Connect function for the database. This wraps around the
        % 'database' function in Matlab to connect to a PostgreSQL
        % database. Connection parameters are set in the constructor when
        % instantiating the class.
        function obj = connect(obj)
            
            fprintf(['Connecting to ' obj.database_name ' database...\n'])
            %%url = 'jdbc:postgresql://host:port/dbname';
            url = ['jdbc:postgresql://' obj.ip_address ':' num2str(obj.port) '/' obj.database_name];
            obj.db_connection = database(obj.database_name,obj.username,obj.password,obj.driver,url);
            %%connect to databse 
            %conn = database(databasename,username,password,driver,url);
            if strncmp(obj.db_connection.Message,'JDBC Driver Error: No suitable driver found',42)     % Database connection status message
                error('MyComponent:incorrectDriver','JDBC Driver Error: No suitable driver found! \nPlease run javaclasspath to check if the JDBC java driver path has been added.\n')
            elseif strncmp(obj.db_connection.Message,'JDBC Driver Error: The connection attempt failed.',48)
                error('MyComponent:incorrectNetwork','JDBC Driver Error: The connection attempt failed! \nPlease check your VPN or Internet connection!\n')
            elseif isempty(obj.db_connection.Message)
                fprintf(['Connected to ' obj.database_name ' database!\n'])
            else
                fprintf(['The connection status is ' obj.db_connection.Message ' !\n'])
            end
            
        end
        
        function fields = parseFields(fields)
           
            for i = 1:length(fields)
               
%                 if strcmp(fields{i}, 'latitude')
%                    
%                     fields{i} = 'ST_Y(geography::geometry) AS latitude';
%                     
%                 end
%                 
%                 if strcmp(fields{i}, 'longitude')
%                    
%                     fields{i} = 'ST_X(geography::geometry) AS longitude';
%                     
%                 end
                
            end
            
        end
        
        % find all the fields given a table
        function fields = allFields(obj,table)
            table_info = sqlfind(obj.db_connection,table);
            fields = table_info(strcmp(table_info.Type,'TABLE'),:).Columns;
            fields = fields{1};
        end
        
        
        % =============================== function convertFromTableToStruct===============================
        % purpose:      convert From Table To Struct with field of array(numeric, cell)
        % intput:       result = data of table format, format: table
        %               column_names = columns of the table
        % output:       result_struct = struct format with field of array(numeric, cell)   
        function result_struct = convertFromTableToStruct(column_names, result)
            
            for i = 1:length(column_names)
%                 result.(column_names{i}) = table2array(result.(column_names{i}));
                if iscell(result.(column_names{i}))
                    result_struct.(column_names{i}) = result.(column_names{i}); %cell2mat
                else % isnumeric
                    result_struct.(column_names{i}) = result.(column_names{i});
                end
            end
            
        end
        
    end %%end of method(static)
    
end