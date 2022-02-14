function [System_Data_Nodes, System_Data_Lines]=renumerate(System_Data_Nodes, System_Data_Lines)

    nb = length(System_Data_Nodes(:,1));
    renum = [System_Data_Nodes(:,1),zeros(length(System_Data_Nodes(:,1)),1)];
    slack = find(System_Data_Nodes(:,2) == 1);
    renum(slack,2) = 1;
    renum(find(System_Data_Nodes(:,2) ~= 1),2) = 2:nb;
    
    System_Data_Nodes(:,1) = renum(:,2);
    
    for i = 1:length(System_Data_Nodes(:,1))
        node = renum(i,1);
        new_node = renum(i,2);
        [r,c] = find(System_Data_Lines(:,1:2) == node);
        for j = 1:length(r)
            System_Data_Lines(r(j),c(j)) = new_node;
        end
    end

end