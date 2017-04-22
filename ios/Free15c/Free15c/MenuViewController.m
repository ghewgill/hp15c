//
//  MenuViewController.m
//  Free15c
//
//  Created by Greg Hewgill on 28/11/12.
//  Copyright (c) 2012 Greg Hewgill. All rights reserved.
//

#import "MenuViewController.h"

@interface MenuViewController ()

@end

@implementation MenuViewController {
    UIWebView *core;
    UISegmentedControl *decimal_style;
}

- (id)initWithCore:(UIWebView *)core_
{
    self = [super init];
    if (self) {
        core = core_;
    }
    return self;
}

- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil
{
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Custom initialization
    }
    return self;
}

- (void)loadView
{
    UIView *frame = [UIView new];
    frame.backgroundColor = [UIColor lightGrayColor];
    
    UINavigationBar *navbar = [UINavigationBar new];
    navbar.translatesAutoresizingMaskIntoConstraints = NO;
    [frame addSubview:navbar];
    
    UINavigationItem *item = [[UINavigationItem alloc] initWithTitle:@"Free 15C"];
    item.rightBarButtonItem = [[UIBarButtonItem alloc] initWithTitle:@"Done" style:UIBarButtonItemStyleDone target:self action:@selector(done)];
    [navbar pushNavigationItem:item animated:NO];

    UILabel *decimal_label = [UILabel new];
    decimal_label.translatesAutoresizingMaskIntoConstraints = NO;
    decimal_label.text = @"Decimal style";
    decimal_label.backgroundColor = [UIColor lightGrayColor];
    decimal_label.font = [UIFont boldSystemFontOfSize:16];
    [frame addSubview:decimal_label];
    
    decimal_style = [[UISegmentedControl alloc] initWithItems:@[@"1,234.56", @"1.234,56"]];
    decimal_style.translatesAutoresizingMaskIntoConstraints = NO;
    decimal_style.selectedSegmentIndex = 0;
    [decimal_style addTarget:nil action:@selector(decimalStyleChanged) forControlEvents:UIControlEventValueChanged];
    [frame addSubview:decimal_style];

    UIWebView *about = [UIWebView new];
    about.translatesAutoresizingMaskIntoConstraints = NO;
    about.backgroundColor = [UIColor lightGrayColor];
    [about loadHTMLString:[NSString stringWithContentsOfFile:[[NSBundle mainBundle] pathForResource:@"about.html" ofType:nil] encoding:NSUTF8StringEncoding error:nil] baseURL:nil];
    about.hidden = YES;
    about.delegate = self;
    [frame addSubview:about];

    [frame addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"V:|[navbar(32)]-20-[decimal_style(44)]-[about]|" options:NSLayoutFormatDirectionLeadingToTrailing metrics:nil views:NSDictionaryOfVariableBindings(navbar, decimal_style, about)]];

    [frame addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[navbar]|" options:NSLayoutFormatDirectionLeadingToTrailing metrics:nil views:NSDictionaryOfVariableBindings(navbar)]];

    [frame addConstraint:[NSLayoutConstraint constraintWithItem:decimal_label attribute:NSLayoutAttributeCenterY relatedBy:NSLayoutRelationEqual toItem:decimal_style attribute:NSLayoutAttributeCenterY multiplier:1 constant:0]];

    [frame addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|-[decimal_label]-[decimal_style]" options:NSLayoutFormatDirectionLeadingToTrailing metrics:nil views:NSDictionaryOfVariableBindings(decimal_label, decimal_style)]];

    [frame addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[about]|" options:NSLayoutFormatDirectionLeadingToTrailing metrics:nil views:NSDictionaryOfVariableBindings(about)]];
    
    self.view = frame;
}

- (void)viewDidLoad
{
    [super viewDidLoad];
	// Do any additional setup after loading the view.
    [UIApplication sharedApplication].statusBarHidden = YES;
    if ([[core stringByEvaluatingJavaScriptFromString:@"DecimalSwap"] isEqualToString:@"1"]) {
        decimal_style.selectedSegmentIndex = 1;
    }
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation
{
    return interfaceOrientation == UIInterfaceOrientationLandscapeLeft
        || interfaceOrientation == UIInterfaceOrientationLandscapeRight;
}

- (BOOL)webView:(UIWebView *)webView shouldStartLoadWithRequest:(NSURLRequest *)request navigationType:(UIWebViewNavigationType)navigationType
{
    if ([request.URL.scheme isEqualToString:@"http"]) {
        [[UIApplication sharedApplication] openURL:request.URL];
        return NO;
    }
    return YES;
}

- (void)webViewDidFinishLoad:(UIWebView *)webView
{
    webView.hidden = NO;
}

- (void)decimalStyleChanged
{
    [core stringByEvaluatingJavaScriptFromString:[NSString stringWithFormat:@"DecimalSwap = %d; update_display();", (int)decimal_style.selectedSegmentIndex]];
    [[NSUserDefaults standardUserDefaults] setInteger:decimal_style.selectedSegmentIndex forKey:@"DecimalStyle"];
    [[NSUserDefaults standardUserDefaults] synchronize];
}

- (void)done
{
    [UIApplication sharedApplication].statusBarHidden = NO;
    [self dismissViewControllerAnimated:YES completion:nil];
}

@end
