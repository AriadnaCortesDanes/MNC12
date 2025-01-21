function interp_gadget(callbak,~)
% Explore the Runge phenomenon -- does polynomial interpolation converge?
% The f popup chooses among three functions:
%    f(x) = 1/(1+c*x^2),
%    f(x) = exp(-c*x*2),
%    f(x) = abs(c*x).
% The + and - buttons increment or decrement three parameters:
%    c: coefficient of x^2,
%    n: number of interpolation points,
%    w: weighted average of equally spaced and Chebyshev points.
% As n increases, does the polynomial interpolant converge to f(x)?
%
% References:
%    Cleve's Corner: https://blogs.mathworks.com/cleve/2018/12/10
%    Wikipedia: https://en.wikipedia.org/wiki/Runge's_phenomenon
%   Copyright 2018 Cleve Moler
%   Copyright 2018 The MathWorks, Inc.
    if nargin == 0 || isequal(get(callbak,'string'),'reset')
        
        % Initialize
        
        c = 25;   % Coefficient in 1/(1+c*x^2) or exp(-c*x^2).
        n = 9;    % Number of interpolation points.
        w = 0.0;  % Weight between equal and Chebyshev spacing.
        
        x = -1:2/(n-1):1;  % Interpolation points.
        u = -1:1/128:1;    % Plot sampling points.     
        
        p = plots(x,u);    
        b = buttons;
        q = labels(c,n,w);
        callbak = false;
        set(gcf,'userdata',{c,n,w,p,q,b,u})
        
    else
             
        % Callbacks
        
        ud = get(gcf,'userdata');
        [c,n,w,p,q,b,u] = deal(ud{:});
        
    end
    
    % Update parameters
    
    [c,n,w] = update(c,n,w,b,q,callbak);
    
    % Which function?
    
    switch get(b(7),'value')
        case {1,2}
            f = @(x) 1./(1+c*x.^2);
        case 3
            f = @(x) exp(-c*x.^2);
        case 4
            f = @(x) abs(c*x);
    end
    
    % Interpolation points
    
    xeq = -1:2/(n-1):1;
    xch = cos((n-1/2:-1:1/2)/n*pi);
    x = (1-w)*xeq + w*xch;
    
    % Polynomial interpolation
    % Interpolate f(x), sampled at x, evaluated at u.
    
    y = polyinterp(x,f(x),u);
    
    % Update plots
    
    set(p(1),'xdata',x,'ydata',f(x))  % Blue circles.
    set(p(2),'xdata',u,'ydata',f(u))  % Thin green line.
    set(p(3),'xdata',u,'ydata',y)     % Interpolant, heavy blue line.
    update_caption(p(4),w,c)
    update_title(p(5),f,c)
    
    % Save parameters and handles
    
    set(gcf,'userdata',{c,n,w,p,q,b,u})
    
% -------------------------------------------------------------
    function p = plots(x,u)
        % Initialize three plots, caption and title.
        clf
        shg
        pos = [.1300 .2100 .7750 .7150];
        ax = axes('position',pos);
        p = plot(x,0*x,'o', ...
                 u,0*u,'-', ...
                 u,0*u,'-');
        set(p(1),'color',[0 0 0.75],'linewidth',2)  % Dark blue
        set(p(2),'color',[0 0.5 0])   % Dark green
        set(p(3),'color',get(p(1),'color'),'linewidth',2)
        set(ax,'xlim',[-1.2 1.2],'ylim',[-.4 1.4])
        p(4) = text(0.0,-0.3,'','horiz','center');  % Caption
        p(5) = title('');
    end % plots
    function b = buttons
        % Initialize six +/- pushbuttons, 'f' popup, ...
        % and info and reset buttons.
        
        b = zeros(1,6);
        for k = 1:6 
            if mod(k,2) == 1
                xp = .13*k+.02;
                s = '-';
            else
                xp = .13*k;
                s = '+';
            end
            b(k) = uicontrol('style','pushbutton', ...
                'units','normalized', ...
                'position',[xp,.02,.07,.07], ...
                'string',s, ... 
                'fontsize',14, ...
                'callback',@interp_gadget);
        end
        
        % Select function
        flist = {'f','1/(1+c*x^2','exp(-c*x^2)','abs(c*x)'};
        b(7) = uicontrol('style','popup', ...
            'units','normalized', ...
            'position',[.92,.50,.07,.07], ...
            'string',flist, ...
            'fontsize',10, ...
            'fontweight','bold', ...
            'callback',@interp_gadget); 
        
        % info button
        
        infocb = @(~,~) ...
             web('http://blogs.mathworks.com/cleve/2018/12/10');
        uicontrol('style','pushbutton', ...
            'units','normalized', ...
            'position',[.03,.02,.07,.07], ...
            'string','info', ... 
            'fontsize',10, ...
            'fontweight','bold', ...
            'callback',infocb);
        
        % reset button
        
        uicontrol('style','pushbutton', ...
            'units','normalized', ...
            'position',[.91,.02,.07,.07], ...
            'string','reset', ... 
            'fontsize',10, ...
            'fontweight','bold', ...
            'callback',@interp_gadget);
    end % buttons
    function q = labels(c,n,w)
        % Initialize three labels.
        q = zeros(1,3);
        for k = 1:3
            if k == 1
                s = sprintf('c = %3.1f',c);
            elseif k == 2
                s = sprintf('n = %2d',n);
            else
                s = sprintf('w = %2.1f',w);
            end
            q(k) = uicontrol('style','text', ...
                'units','normalized', ...
                'position',[.26*k-.08,.10,.14,.05], ...
                'string',s, ...
                'fontsize',12, ...
                'fontweight','bold');    
        end     
    end % labels
 
    function [c,n,w] = update(c,n,w,b,q,callbak)
        % Respond to callback.
        switch callbak
            case b(1)
                c = c-delc(b);
            case b(2)
                c = c+delc(b);
            case b(3)
                n = n-1;
            case b(4)
                n = n+1;
            case b(5)
                w = w-0.1;
            case b(6)
                w = w+0.1;
            case b(7)
                if get(b(7),'value') < 4
                    c = 25;
                else
                    c = 1;
                end
        end
        set(q(1),'string',sprintf('c = %3.1f',c));
        set(q(2),'string',sprintf('n = %2.0f',n));
        set(q(3),'string',sprintf('w = %2.1f',w));
    end % update
    function del = delc(b)
        if get(b(7),'value') < 4
            del = 1;
        else 
            del = .1;
        end
    end % delc
    function update_caption(caption,w,c)
        if abs(w) < eps
            s = 'equally spaced';
        elseif abs(1-w) < eps
            s = 'chebyshev';
        elseif w < 1
            s = sprintf('%2.1f*cheby + %2.1f*equal',w,1-w);
        else
            s = sprintf('%2.1f*cheby - %2.1f*equal',w,w-1);
        end
        if c > 0
            ypos = -0.3;
            set(gca,'ylim',[-0.4 1.4])
        else
            ypos = -7;
            set(gca,'ylim',[-8 10])
        end
        set(caption,'string',s,'position',[0.0 ypos 0.0])
    end % update_caption
    function update_title(titl,f,c)
        F = func2str(f);
        F = F(5:end);  % Remove leading '(x)'.
        F = regexprep(F,'\.',''); % Remove dots.
        if c == 0
            s = '1';
        elseif c == 1
            s = regexprep(F,'c\*','');
        elseif c > 0
            if c >= 10
                form = regexprep(F,'c','%2.0f');
            else
                form = regexprep(F,'c','%3.1f');
            end
            s = sprintf(form,c);
        elseif c == -1 && any(F=='-')
            s = regexprep(F,'-c\*','');
        elseif c == -1   
            s = regexprep(F,'+c\*','-');
        elseif any(F=='-')
            form = regexprep(F,'-c','%2.0f');
            s = sprintf(form,-c);
        else
            if c <= -10
                form = regexprep(F,'c','%2.0f');
            else
                form = regexprep(F,'c','%3.1f');
            end
            s = sprintf(form,-c);
        end
        set(titl,'string',s)   
    end % update_title
    function v = polyinterp(x,y,u)
        %POLYINTERP  Polynomial interpolation.
        %   v = POLYINTERP(x,y,u) computes v(j) = P(u(j)) where P is the
        %   polynomial of degree d = length(x)-1 with P(x(i)) = y(i).
        % Use Lagrangian representation.
        % Evaluate at all elements of u simultaneously.
        n = length(x);
        v = zeros(size(u));
        for k = 1:n
           z = ones(size(u));
           for j = [1:k-1 k+1:n]
              z = (u-x(j))./(x(k)-x(j)).*z;
           end
           v = v + z*y(k);
        end
    end % polyinterp
end % interp_gadget