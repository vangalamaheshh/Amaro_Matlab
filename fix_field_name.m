function fn=fix_field_name(str)

fn=str;
fn(fn==' ')='_';
fn(fn=='(')='_';
fn(fn==')')='_';
fn(fn=='-')='_';
